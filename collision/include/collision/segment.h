/**
 * @file segment.h
 * @author khalilhenoud@gmail.com
 * @brief
 * @version 0.1
 * @date 2023-06-10
 *
 * @copyright Copyright (c) 2023
 *
 */
#ifndef SEGMENT_COLLISION_H
#define SEGMENT_COLLISION_H

#ifdef __cplusplus
extern "C" {
#endif

#include <collision/internal/module.h>
#include <collision/segment_classification.h>
#include <math/c/segment.h>


// NOTE: this function guarantees that the first point in 'out' belongs to
// 'first' and vice versa (full overlap cases are the exception for obvious
// reasons).
COLLISION_API
segments_classification_t
classify_segments(
  const segment_t *first,
  const segment_t *second,
  segment_t *out);

#ifdef __cplusplus
}
#endif

#endif