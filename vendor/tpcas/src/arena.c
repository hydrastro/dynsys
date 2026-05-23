#include "arena.h"
#include "lib/allocators.h"

#include <stdint.h>
#include <string.h>

struct arena_chunk {
    ds_arena_t    storage;
    arena_chunk_t *next;
    size_t         capacity;
};

static uintptr_t align_up(uintptr_t value, size_t align) {
    if (align <= 1) return value;
    uintptr_t mask = (uintptr_t)align - 1u;
    return (value + mask) & ~mask;
}

static arena_chunk_t *arena_new_chunk(arena_t *a, size_t capacity) {
    arena_chunk_t *chunk = ds_context_calloc(&a->backing, 1, sizeof *chunk);
    if (!chunk) return NULL;
    if (ds_arena_create(&chunk->storage, capacity, &a->backing) != DS_OK) {
        ds_context_free(&a->backing, chunk);
        return NULL;
    }
    chunk->capacity = capacity;
    return chunk;
}

static void arena_free_chunk(arena_t *a, arena_chunk_t *chunk) {
    if (!chunk) return;
    ds_arena_destroy(&chunk->storage, &a->backing);
    ds_context_free(&a->backing, chunk);
}

static void *chunk_alloc_align(arena_chunk_t *chunk, size_t bytes, size_t align) {
    ds_context_t alloc;
    unsigned char *raw;
    size_t request = bytes;

    if (bytes == 0) bytes = request = 1;
    if (align > sizeof(void *)) {
        if (bytes > (size_t)-1 - align) return NULL;
        request = bytes + align;
    }

    ds_context_init(&alloc);
    ds_context_use_arena(&alloc, &chunk->storage);
    raw = ds_context_calloc(&alloc, 1, request);
    if (!raw) return NULL;
    if (align <= sizeof(void *)) return raw;
    return (void *)align_up((uintptr_t)raw, align);
}

void arena_init(arena_t *a, size_t chunk_size) {
    if (!a) return;
    ds_context_init(&a->backing);
    a->chunk_size = chunk_size ? chunk_size : 65536;
    a->head = arena_new_chunk(a, a->chunk_size);
}

void *arena_alloc_align(arena_t *a, size_t bytes, size_t align) {
    if (!a) return NULL;
    if (align == 0) align = 1;
    if (!a->head) {
        a->head = arena_new_chunk(a, a->chunk_size ? a->chunk_size : 65536);
        if (!a->head) return NULL;
    }

    void *p = chunk_alloc_align(a->head, bytes, align);
    if (p) return p;

    size_t need = bytes + align;
    size_t cap = a->chunk_size > need ? a->chunk_size : need * 2;
    arena_chunk_t *chunk = arena_new_chunk(a, cap);
    if (!chunk) return NULL;
    chunk->next = a->head;
    a->head = chunk;
    return chunk_alloc_align(a->head, bytes, align);
}

void *arena_alloc(arena_t *a, size_t bytes) {
    return arena_alloc_align(a, bytes, sizeof(void *));
}

char *arena_strndup(arena_t *a, const char *s, size_t n) {
    if (!s) return NULL;
    char *p = arena_alloc_align(a, n + 1, 1);
    if (!p) return NULL;
    memcpy(p, s, n);
    p[n] = '\0';
    return p;
}

char *arena_strdup(arena_t *a, const char *s) {
    return s ? arena_strndup(a, s, strlen(s)) : NULL;
}

void arena_reset(arena_t *a) {
    if (!a || !a->head) return;
    arena_chunk_t *chunk = a->head->next;
    while (chunk) {
        arena_chunk_t *next = chunk->next;
        arena_free_chunk(a, chunk);
        chunk = next;
    }
    a->head->next = NULL;
    ds_arena_reset(&a->head->storage);
}

void arena_destroy(arena_t *a) {
    if (!a) return;
    arena_chunk_t *chunk = a->head;
    while (chunk) {
        arena_chunk_t *next = chunk->next;
        arena_free_chunk(a, chunk);
        chunk = next;
    }
    a->head = NULL;
    a->chunk_size = 0;
}

size_t arena_used(const arena_t *a) {
    size_t total = 0;
    for (arena_chunk_t *chunk = a ? a->head : NULL; chunk; chunk = chunk->next) {
        total += ds_arena_used(&chunk->storage);
    }
    return total;
}

size_t arena_remaining(const arena_t *a) {
    size_t total = 0;
    for (arena_chunk_t *chunk = a ? a->head : NULL; chunk; chunk = chunk->next) {
        total += ds_arena_remaining(&chunk->storage);
    }
    return total;
}
