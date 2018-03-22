/***********************************************************
        

        Starter code for Assignment 3

        Implements light_source.h

***********************************************************/

#include "light_source.h"
#include <cmath>

void PointLight::shade(Ray3D &ray) {
    // TODO: implement this function to fill in values for ray.col
    // using phong shading.  Make sure your vectors are normalized, and
    // clamp colour values to 1.0.
    //
    // It is assumed at this point that the intersection information in ray
    // is available.  So be sure that traverseScene() is called on the ray
    // before this function.

    // Getting values from ray
    Intersection intersection = ray.intersection;
    Material* mat = intersection.mat;

    // normalized normal at intersection point
    Vector3D normal = intersection.normal;
    normal.normalize();

    // normalized light ray
    Vector3D light_ray = get_position() - intersection.point;
    light_ray.normalize();

    // direction vector
    Vector3D dir = ray.dir;

    // reflection of light_ray on the normal
    Vector3D R = 2.0 * (normal.dot(light_ray)) * normal - light_ray;

    Color ambient_term = mat->ambient * col_ambient;
    Color diffuse_term = fmax(light_ray.dot(normal),0) * (mat->diffuse * col_diffuse);
    Color specular_term = pow(fmax(-R.dot(dir), 0), mat->specular_exp) * (mat->specular * col_specular);

    Color col = ray.col;
    col = col + ambient_term + diffuse_term + specular_term;
    col.clamp();

    ray.col = col;
}
