#ifndef ARENA_H
#define ARENA_H

#include <stddef.h>
#include "lib/context.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TPCAS keeps a small project-local arena facade, but the storage policy now
 * comes from the custom data-structures library. Internally each growth chunk
 * is a ds_arena_t, preserving the old “grow by chunks” behavior while reusing
 * the DS allocator/context layer. */

typedef struct arena_chunk arena_chunk_t;

typedef struct {
    arena_chunk_t *head;
    size_t         chunk_size;
    ds_context_t   backing;
} arena_t;

void  arena_init(arena_t *a, size_t chunk_size);
void *arena_alloc(arena_t *a, size_t bytes);
void *arena_alloc_align(arena_t *a, size_t bytes, size_t align);
char *arena_strdup(arena_t *a, const char *s);
char *arena_strndup(arena_t *a, const char *s, size_t n);
void  arena_reset(arena_t *a);
void  arena_destroy(arena_t *a);
size_t arena_used(const arena_t *a);
size_t arena_remaining(const arena_t *a);


#ifdef __cplusplus
}
#endif
#endif
