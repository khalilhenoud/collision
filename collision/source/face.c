/**
 * @file face.c
 * @author khalilhenoud@gmail.com
 * @brief 
 * @version 0.1
 * @date 2023-06-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <math.h>
#include <assert.h>
#include <collision/face.h>
#include <collision/segment.h>
#include <collision/sphere.h>
#include <collision/capsule.h>


void
get_faces_normals(
  const face_t *faces, 
  const uint32_t count, 
  vector3f *normals)
{
  vector3f v1, v2;

  assert(normals != NULL);

  for (uint32_t i = 0; i < count; ++i) {
    vector3f_set_diff_v3f(&v1, &faces[i].points[0], &faces[i].points[1]);
    vector3f_set_diff_v3f(&v2, &faces[i].points[0], &faces[i].points[2]);
    normals[i] = cross_product_v3f(&v1, &v2);
    normalize_set_v3f(normals + i);
  }
}

float
get_point_distance(
  const face_t *face, 
  const vector3f *normal, 
  const point3f *point)
{
  vector3f to_point = diff_v3f(&face->points[0], point);
  return dot_product_v3f(normal, &to_point);
}

point3f
get_point_projection(
  const face_t *face, 
  const vector3f *normal, 
  const point3f *point,
  float *distance)
{
  assert(distance != NULL);
  *distance = get_point_distance(face, normal, point);
  
  {
    vector3f scaled_normal = mult_v3f(normal, *distance);
    point3f projected = diff_v3f(&scaled_normal, point);
    return projected;
  }
}

segment_plane_classification_t
classify_segment_face(
  const face_t *face, 
  const vector3f *normal, 
  const segment_t *segment,
  point3f *intersection,
  float *t)
{
  assert(intersection != NULL);
  assert(t != NULL);

  vector3f delta = diff_v3f(segment->points + 0, segment->points + 1);
  const point3f *A = segment->points + 0;

  {
    float dotproduct = dot_product_v3f(&delta, normal);
    if (IS_ZERO_LP(dotproduct)) {
      // colinear or parallel.
      float distance = get_point_distance(face, normal, A);
      return IS_ZERO_LP(distance) ? 
        SEGMENT_PLANE_COPLANAR : 
        SEGMENT_PLANE_PARALLEL;
    }
  }

  {
    float D = 
      - normal->data[0] * face->points[0].data[0] 
      - normal->data[1] * face->points[0].data[1] 
      - normal->data[2] * face->points[0].data[2];
    
    float denom = 
      normal->data[0] * delta.data[0] + 
      normal->data[1] * delta.data[1] + 
      normal->data[2] * delta.data[2];
    assert(!IS_ZERO_LP(denom));
    *t = (
      - normal->data[0] * A->data[0] 
      - normal->data[1] * A->data[1] 
      - normal->data[2] * A->data[2] - D) / denom;

    // scale the delta by t, and add it to A to get the intersection.
    mult_set_v3f(&delta, *t);
    *intersection = add_v3f(A, &delta);
    
    {
      // nudge to plane, basically error correction.
      float distance = get_point_distance(face, normal, intersection);
      vector3f scaled_normal = *normal;
      scaled_normal = mult_v3f(&scaled_normal, -distance);
      add_set_v3f(intersection, &scaled_normal);
    }

    return (*t <= 1.f && *t >= 0.f) ? 
      SEGMENT_PLANE_INTERSECT_ON_SEGMENT : 
      SEGMENT_PLANE_INTERSECT_OFF_SEGMENT;
  }

  return SEGMENT_PLANE_PARALLEL;
}

coplanar_point_classification_t
classify_coplanar_point_face(
  const face_t *face, 
  const vector3f *normal,
  const point3f *point,
  point3f *closest)
{
  assert(closest != NULL);

  // NOTE: The function will not check if the point is coplanar or not. It is up
  // to the caller to make sure it is. The reason for this is the epsilon value
  // choice, which is unreliable and error prone.

  {
    vector3f v01, v12, v20;
    vector3f v0c, v1c, v2c;
    vector3f c[3];
    float dot[3] = { 0 };

    vector3f_set_diff_v3f(&v01, &face->points[0], &face->points[1]);
    normalize_set_v3f(&v01);
    vector3f_set_diff_v3f(&v0c, &face->points[0], point);
    normalize_set_v3f(&v0c);
    c[0] = cross_product_v3f(&v01, &v0c);
    dot[0] = dot_product_v3f(c + 0, normal);
    
    vector3f_set_diff_v3f(&v12, &face->points[1], &face->points[2]);
    normalize_set_v3f(&v12);
    vector3f_set_diff_v3f(&v1c, &face->points[1], point);
    normalize_set_v3f(&v1c);
    c[1] = cross_product_v3f(&v12, &v1c);
    dot[1] = dot_product_v3f(c + 1, normal);
    
    vector3f_set_diff_v3f(&v20, &face->points[2], &face->points[0]);
    normalize_set_v3f(&v20);
    vector3f_set_diff_v3f(&v2c, &face->points[2], point);
    normalize_set_v3f(&v2c);
    c[2] = cross_product_v3f(&v20, &v2c);
    dot[2] = dot_product_v3f(c + 2, normal);

    if (
      (
        dot[0] <= EPSILON_FLOAT_LOW_PRECISION && 
        dot[1] <= EPSILON_FLOAT_LOW_PRECISION && 
        dot[2] <= EPSILON_FLOAT_LOW_PRECISION) ||
      (
        dot[0] >= -EPSILON_FLOAT_LOW_PRECISION && 
        dot[1] >= -EPSILON_FLOAT_LOW_PRECISION && 
        dot[2] >= -EPSILON_FLOAT_LOW_PRECISION)) {
      *closest = *point;
      return COPLANAR_POINT_ON_OR_INSIDE;
    } else {
      vector3f vec, p[3];
      float dist[3];
      p[0] = closest_point_on_segment_loose(
        point, 
        &face->points[0], 
        &face->points[1]);
      p[1] = closest_point_on_segment_loose(
        point, 
        &face->points[1], 
        &face->points[2]);
      p[2] = closest_point_on_segment_loose(
        point, 
        &face->points[2], 
        &face->points[0]);
      vector3f_set_diff_v3f(&vec, p + 0, point);
      dist[0] = length_v3f(&vec);
      vector3f_set_diff_v3f(&vec, p + 1, point);
      dist[1] = length_v3f(&vec);
      vector3f_set_diff_v3f(&vec, p + 2, point);
      dist[2] = length_v3f(&vec);

      {
        float minimum = dist[0];
        uint32_t index = 0;

        for (uint32_t i = 1; i < 3; ++i) {
          if (dist[i] < minimum) {
            minimum = dist[i];
            index = i;
          }
        }

        *closest = *(p + index);
      }

      return COPLANAR_POINT_OUTSIDE;
    }
  }
}

point_halfspace_classification_t
classify_point_halfspace(
  const face_t *face, 
  const vector3f *normal, 
  const point3f *point)
{
  float distance = get_point_distance(face, normal, point);
  if (IS_ZERO_LP(distance)) 
    return POINT_ON_PLANE;
  else if (distance > 0.f)
    return POINT_IN_POSITIVE_HALFSPACE;
  else 
    return POINT_IN_NEGATIVE_HALFSPACE;
}

face_t
get_extended_face(
  const face_t* face,
  float radius)
{
  assert(face);

  {
    face_t augmented;
    float angles[3];
    float dot_product;
    float k;
    vector3f offset;
    vector3f avec[2];
    uint32_t vec0[3][2] = { {0, 1}, {1, 2}, {2, 0} };
    uint32_t vec1[3][2] = { {0, 2}, {1, 0}, {2, 1} };
    for (uint32_t i = 0; i < 3; ++i) {
      vector3f_set_diff_v3f(
        avec + 0, &face->points[vec0[i][0]], &face->points[vec0[i][1]]);
      normalize_set_v3f(avec + 0);
      vector3f_set_diff_v3f(
        avec + 1, &face->points[vec1[i][0]], &face->points[vec1[i][1]]);
      normalize_set_v3f(avec + 1);

      dot_product = dot_product_v3f(avec + 0, avec + 1);
      angles[i] = acosf(fabs(dot_product));

      k = radius / sinf(angles[i]);
      vector3f_copy(augmented.points + i, face->points + i);
      mult_set_v3f(avec + 0, -k);
      mult_set_v3f(avec + 1, -k);
      vector3f_copy(&offset, avec + 0);
      add_set_v3f(&offset, avec + 1);
      add_set_v3f(augmented.points + i, &offset); 
    }

    return augmented;
  }
}

static
void
calculate_sphere_faceplane_penetration(
  const sphere_t* sphere,
  const face_t* face,
  const vector3f* normal,
  const float distance,
  vector3f* penetration)
{
  float scale = sphere->radius - distance;
  *penetration = *normal;
  mult_set_v3f(penetration, scale);
}

static
void
calculate_sphere_face_penetration(
  const sphere_t* sphere,
  const face_t* face,
  const vector3f* normal,
  const float distance,
  vector3f* closest,
  vector3f* penetration)
{
  vector3f from_to = diff_v3f(closest, &sphere->center);
  float scale = sphere->radius - length_v3f(&from_to);
  normalize_set_v3f(&from_to);
  mult_set_v3f(&from_to, scale);
  *penetration = from_to;
}

sphere_face_classification_t
classify_sphere_face(
  const sphere_t* sphere,
  const face_t* face,
  const vector3f* normal,
  const int32_t as_plane,
  vector3f* penetration)
{
  float distance = 0.f;
  vector3f projected = get_point_projection(
    face, 
    normal, 
    &sphere->center, 
    &distance);
  
  if (fabs(distance) > sphere->radius)
    return SPHERE_FACE_NO_COLLISION;

  {
    // check if the projected sphere center is in the face.
    vector3f closest_on_face;
    coplanar_point_classification_t classification;
    classification = classify_coplanar_point_face(
      face, 
      normal, 
      &projected, 
      &closest_on_face);

    if (classification != COPLANAR_POINT_ON_OR_INSIDE) {
      vector3f from_to = diff_v3f(&closest_on_face, &sphere->center);
      if (length_v3f(&from_to) > sphere->radius)
        return SPHERE_FACE_NO_COLLISION;
    }

    // guaranteed collision at this point.
    if (as_plane) {
      calculate_sphere_faceplane_penetration(
        sphere, face, normal, distance, penetration);
    } else {
      calculate_sphere_face_penetration(
        sphere, 
        face, 
        normal, 
        distance, 
        classification == COPLANAR_POINT_ON_OR_INSIDE ? 
          &projected : &closest_on_face, 
        penetration);
    }

    return classification == COPLANAR_POINT_ON_OR_INSIDE ? 
      (IS_ZERO_LP(distance) ? 
        SPHERE_FACE_COLLIDES_SPHERE_CENTER_ON_FACE : 
        SPHERE_FACE_COLLIDES) : 
        SPHERE_FACE_COLLIDES;
  }
}

static
void
extrude_capsule_along_face_normal(
  const capsule_t* capsule,
  const face_t* face,
  const vector3f* normal,
  const segment_t* segment,
  const point3f* intersection,
  vector3f* penetration)
{
  // shift the segment points by -radius along the normal.
  int32_t index = -1;
  segment_t shifted;
  vector3f to_shift = *normal;
  mult_set_v3f(&to_shift, -capsule->radius);
  shifted.points[0] = add_v3f(segment->points + 0, &to_shift);
  shifted.points[1] = add_v3f(segment->points + 1, &to_shift);

  // find the point with the most negative distance from the face.
  {
    float distance[] = { 0.f, 0.f };
    distance[0] = get_point_distance(face, normal, shifted.points + 0);
    distance[1] = get_point_distance(face, normal, shifted.points + 1);
    index = (distance[0] < distance[1]) ? 0 : 1;

    assert(index != -1);
    assert(
      distance[index] < EPSILON_FLOAT_LOW_PRECISION && 
      "The shifted point has to be in the negative space!");
  }

  {
    // find the projection of that point on the face.
    float distance = 0.f;
    point3f closest;
    point3f projected = get_point_projection(
      face, normal, shifted.points + index, &distance);
    coplanar_point_classification_t classification = 
      classify_coplanar_point_face(face, normal, &projected, &closest);
    
    if (classification == COPLANAR_POINT_ON_OR_INSIDE) {
      distance = fabs(distance);
      *penetration = *normal;
      mult_set_v3f(penetration, distance);
    } else {
      point3f closest_on_segment;
      segment_t to_intersect_face;
      to_intersect_face.points[0] = *intersection;
      to_intersect_face.points[1] = projected;

      {
        segment_t face_segments[3] = { 
          {face->points[0], face->points[1]}, 
          {face->points[1], face->points[2]}, 
          {face->points[2], face->points[0]}};
        segment_t result;
        segments_classification_t classify;
        int32_t found = 0;

        for (int32_t i = 0; i < 3; ++i) {
          // TODO: write a separate classify_coplanar_segments function and use
          // it here. This is to avoid imprecision errors. Then assert if found
          // is 0. You can remove the if test after this one.
          classify = classify_segments(
            &to_intersect_face, face_segments + i, &result);
          if (classify == SEGMENTS_COPLANAR_INTERSECT) {
            closest_on_segment = result.points[0];
            found = 1;
            break;
          }
        }

        if (found == 0) {
          // NOTE: this does happen when the closest point on the triangle is
          // very close to the projected point. The segment in that case is very
          // close to perpendicular to the plane. This manifest by the
          // classify_segments function indicating the segments are not coplanar
          // by very small distnaces. The easiest way to deal with this is to
          // write a separate fucntion that assumes the segments are coplanar.
          distance = fabs(distance);
          *penetration = *normal;
          mult_set_v3f(penetration, distance);
          return;
        }

        {
          vector3f closest_to_intersection = diff_v3f(
            intersection, &closest_on_segment);
          float length = length_v3f(&closest_to_intersection);
          *penetration = *normal;
          mult_set_v3f(penetration, length);
        }
      }
    }
  }
}

capsule_face_classification_t
classify_capsule_face(
  const capsule_t* capsule,
  const face_t* face,
  const vector3f* normal,
  const int32_t as_plane,
  vector3f* penetration,
  point3f* sphere_center)
{
  segment_t segment;
  get_capsule_segment(capsule, &segment);
  
  assert(penetration != NULL);

  {
    // get the intersection of the segment with the plane.
    float t;
    point3f intersection;
    segment_plane_classification_t classify_segment = 
      classify_segment_face(face, normal, &segment, &intersection, &t);
    
    if (classify_segment == SEGMENT_PLANE_COPLANAR) {
      coplanar_point_classification_t p_classify[2];
      point3f closest[2];
      p_classify[0] = classify_coplanar_point_face(
        face, normal, segment.points + 0, closest + 0);
      p_classify[1] = classify_coplanar_point_face(
        face, normal, segment.points + 1, closest + 1);

      if (
        p_classify[0] == COPLANAR_POINT_ON_OR_INSIDE || 
        p_classify[1] == COPLANAR_POINT_ON_OR_INSIDE) {
        *penetration = *normal;
        mult_set_v3f(penetration, capsule->radius);
        return CAPSULE_FACE_COLLIDES_CAPSULE_AXIS_COPLANAR_FACE;
      } else {
        {
          uint32_t jval[] = { 1, 2, 0, (uint32_t)-1 };

          for (uint32_t i = 0, j = jval[i]; i < 3; ++i, j = jval[i]) {
            segments_classification_t seg_classify;
            segment_t face_segment, out_val;
            face_segment.points[0] = face->points[i];
            face_segment.points[1] = face->points[j];
            seg_classify = classify_segments(
              &face_segment, &segment, &out_val);

            assert(
              seg_classify != SEGMENTS_IDENTICAL &&
              seg_classify != SEGMENTS_COLINEAR_OVERLAPPING &&
              seg_classify != SEGMENTS_COLINEAR_OVERLAPPING_AT_POINT &&
              seg_classify != SEGMENTS_DISTINCT &&
              "the segments cannot be distinct if already coplanar!");

            if (
              seg_classify == SEGMENTS_COPLANAR_INTERSECT || 
              seg_classify == SEGMENTS_COLINEAR_FULL_OVERLAP) {
              *penetration = *normal;
              mult_set_v3f(penetration, capsule->radius);
              return CAPSULE_FACE_COLLIDES_CAPSULE_AXIS_COPLANAR_FACE;
            } 
          }
        }

        // both points are outside the face, this devolves into sphere-face
        // classify. sphere-center would be the closest point on the segment.
        {
          point3f on_capsule[3];
          uint32_t index = 0;
          on_capsule[0] = closest_point_on_segment(
            face->points + 0, &segment);
          on_capsule[1] = closest_point_on_segment(
            face->points + 1, &segment);
          on_capsule[2] = closest_point_on_segment(
            face->points + 2, &segment);

          {
            vector3f diff = diff_v3f(face->points + 0, on_capsule + 0);
            float min_dis = length_squared_v3f(&diff), temp = 0.f;
            for (uint32_t i = 1; i < 3; ++i) {
              diff = diff_v3f(face->points + i, on_capsule + i);
              temp = length_squared_v3f(&diff);
              if (temp < min_dis) {
                min_dis = temp;
                index = i;
              }
            }
          }

          {
            sphere_t sphere = { on_capsule[index], capsule->radius };
            sphere_face_classification_t classify_sphere = 
              classify_sphere_face(
                &sphere, face, normal, as_plane, penetration);

            *sphere_center = sphere.center;
            return classify_sphere == SPHERE_FACE_NO_COLLISION ?
              CAPSULE_FACE_NO_COLLISION : 
              CAPSULE_FACE_COLLIDES_CAPSULE_AXIS_COPLANAR_FACE;
          }
        }
      }
    } else if (classify_segment == SEGMENT_PLANE_PARALLEL) {
      // classify all 3 face segments and the actual segment. Use the minimum
      // connecting branch.
      point3f on_capsule_axis;
      segments_classification_t seg_classify[3];
      segment_t results[3];
      vector3f vec;
      float dist;
      float min_dis = FLT_MAX;
      uint32_t jval[] = { 1, 2, 0, (uint32_t)-1 };

      for (uint32_t i = 0, j = jval[i]; i < 3; ++i, j = jval[i]) {
        segment_t face_segment;
        face_segment.points[0] = face->points[i];
        face_segment.points[1] = face->points[j];
        seg_classify[i] = classify_segments(
          &face_segment, &segment, results + i);
        vector3f_set_diff_v3f(&vec, results[i].points, results[i].points + 1);
        dist = length_squared_v3f(&vec);
        if (dist < min_dis) {
          min_dis = dist;
          on_capsule_axis = results[i].points[1];
        }
      }

      {
        sphere_t sphere = { on_capsule_axis, capsule->radius };
        sphere_face_classification_t classify_sphere = 
          classify_sphere_face(
            &sphere, face, normal, as_plane, penetration);

        *sphere_center = sphere.center;
        assert(
          classify_sphere != SPHERE_FACE_COLLIDES_SPHERE_CENTER_ON_FACE &&
          "We already know the capsule axis is parallel but not coplanar!!");

        return classify_sphere == SPHERE_FACE_NO_COLLISION ?
          CAPSULE_FACE_NO_COLLISION : CAPSULE_FACE_COLLIDES;
      }
    } else {
      // the segment does intersect the face plane.
      point3f closest;
      coplanar_point_classification_t classify_point = 
        classify_coplanar_point_face(face, normal, &intersection, &closest);

      if (
        classify_segment == SEGMENT_PLANE_INTERSECT_OFF_SEGMENT ||
        classify_point == COPLANAR_POINT_OUTSIDE) {
        // either the intersection point is outside the face or it is inside but
        // outside the boundaries of the segment.
        point3f on_capsule_axis = closest_point_on_segment(&closest, &segment);
        sphere_t sphere = { on_capsule_axis, capsule->radius };
        sphere_face_classification_t classify_sphere = 
          classify_sphere_face(
            &sphere, face, normal, as_plane, penetration);
        *sphere_center = sphere.center;

        if (classify_sphere == SPHERE_FACE_NO_COLLISION)
          return CAPSULE_FACE_NO_COLLISION;
        else {
          if (
            as_plane && 
            classify_segment == SEGMENT_PLANE_INTERSECT_ON_SEGMENT) {
            face_t extended = get_extended_face(face, capsule->radius);
            coplanar_point_classification_t classify_point_extended =
              classify_coplanar_point_face(
                &extended, normal, &intersection, &closest);

            if (classify_point_extended == COPLANAR_POINT_ON_OR_INSIDE) {
              extrude_capsule_along_face_normal(
                capsule, 
                &extended, 
                normal, 
                &segment, 
                &intersection, 
                penetration);
            }
          }

          return CAPSULE_FACE_COLLIDES;
        }
      } else {
        extrude_capsule_along_face_normal(
          capsule, face, normal, &segment, &intersection, penetration);
        return CAPSULE_FACE_COLLIDES_CAPSULE_AXIS_INTERSECTS_FACE;
      }
    }
  }

  return CAPSULE_FACE_NO_COLLISION;
}