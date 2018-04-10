/***********************************************************
        

        Starter code for Assignment 3

        Implementations of functions in raytracer.h,
        and the main function which specifies the scene to be rendered.

***********************************************************/

#include "raytracer.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

const int MAX_REFLECT = 1;
const double ERR = 1e-5;

void Raytracer::traverseScene(Scene &scene, Ray3D &ray) {
    for (size_t i = 0; i < scene.size(); ++i) {
        SceneNode *node = scene[i];

        if (node->obj->intersect(ray, node->worldToModel, node->modelToWorld)) {
            ray.intersection.mat = node->mat;
        }
    }
}

void Raytracer::computeTransforms(Scene &scene) {
    // right now this method might seem redundant. But if you decide to
    // implement scene graph this is where you would propagate transformations
    // to child nodes

    for (size_t i = 0; i < scene.size(); ++i) {
        SceneNode *node = scene[i];

        node->modelToWorld = node->trans;
        node->worldToModel = node->invtrans;
    }
}

void Raytracer::computeShading(Ray3D &ray, LightList &light_list, Scene &scene) {
    for (size_t i = 0; i < light_list.size(); ++i) {
        LightSource *light = light_list[i];

        // light->shade(ray);

        // Each lightSource provides its own shading function.
        // Implement shadows here if needed.
        // Check if there is an object between the light source and the intersection point
        Vector3D light_dir = light->get_position() - ray.intersection.point;
        light_dir.normalize();
        Point3D intersect = ray.intersection.point + (ERR*ray.intersection.normal) + (ERR*light_dir);
        bool blocked = false;

        Ray3D new_ray(intersect, light_dir);
        for (size_t i = 0; i < scene.size(); ++i) {
            SceneNode *node = scene[i];

            if (node->obj->intersect(new_ray, node->worldToModel, node->modelToWorld)) {
                blocked = true;
                break;
            }
        }
        if(blocked){
            ray.col = Color(0.0, 0.0, 0.0);
        } else {
            light->shade(ray);
        }
    }
}

Color Raytracer::shadeRay(Ray3D &ray, Scene &scene, LightList &light_list, int num_reflect) {
    Color col(0.0, 0.0, 0.0);
    traverseScene(scene, ray);

    // Don't bother shading if the ray didn't hit
    // anything.
    if (!ray.intersection.none) {
        computeShading(ray, light_list, scene);
        col = ray.intersection.mat->transparency * ray.col;
    }

    // You'll want to call shadeRay recursively (with a different ray,
    // of course) here to implement reflection/refraction effects.

    // Reflect a ray MAX_REFLECT number of times
    // TODO: Fix implementation 
    if (!ray.intersection.none && ray.intersection.mat->reflective) {
        if(num_reflect < MAX_REFLECT){
            Intersection intersect = ray.intersection;
            Vector3D normal = intersect.normal;
            normal.normalize();

            Vector3D light_ray = -ray.dir;
            Vector3D reflect_dir = 2.0 * normal.dot(light_ray) * normal - light_ray;
            reflect_dir.normalize();
            
            Point3D new_intersect = intersect.point + (ERR * intersect.normal);
            // Shoot a new ray from the intersect point in the direction of reflection
            Ray3D new_ray(new_intersect, reflect_dir);
            
            num_reflect++;
            new_ray.col = shadeRay(new_ray, scene, light_list, num_reflect);
        
            // col = col + shadeRay(new_ray, scene, light_list, num_reflect);
            col = col + (intersect.mat->transparency * intersect.mat->specular * new_ray.col);
        }
    }
    col.clamp();
    return col;
}

void Raytracer::render(Camera &camera, Scene &scene, LightList &light_list,
                       Image &image) {
    computeTransforms(scene);

    Matrix4x4 viewToWorld;
    double factor = (double(image.height) / 2) / tan(camera.fov * M_PI / 360.0);

    viewToWorld = camera.initInvViewMatrix();
    // Construct a ray for each pixel.
    Color color;
    double samples = 5;
    for (int i = 0; i < image.height; i++) {
        for (int j = 0; j < image.width; j++) {
            // Implement anti aliasing with samples*samples number
            // of rays for each pixel
            for (int k = 0; k < samples; k++) {
                for(int l = 0; l < samples; l++) {
                    double rand_width = ((double) rand() / RAND_MAX)/samples + l/samples + j;
                    double rand_height = ((double) rand() / RAND_MAX)/samples + k/samples + i;

                    Point3D origin(0, 0, 0);
                    Point3D imagePlane;

                    imagePlane[0] = (-double(image.width) / 2 + rand_width) / factor;
                    imagePlane[1] = (-double(image.height) / 2 + rand_height) / factor;
                    imagePlane[2] = -1;

                    Ray3D ray;
                    // TODO: Convert ray to world space
                    ray.origin = viewToWorld * origin;
                    ray.dir = viewToWorld * (imagePlane - origin);
                    ray.dir.normalize();
                    
                    // Adding colors together
                    color = color + shadeRay(ray, scene, light_list, 0);
                }
            }
        
            // Sets up ray origin and direction in view space,
            // image plane is at z = -1.
            // Point3D origin(0, 0, 0);
            // Point3D imagePlane;
            // imagePlane[0] = (-double(image.width) / 2 + 0.5 + j) / factor;
            // imagePlane[1] = (-double(image.height) / 2 + 0.5 + i) / factor;
            // imagePlane[2] = -1;

            // Ray3D ray;
            // // TODO: Convert ray to world space
            // ray.origin = viewToWorld * origin;
            // ray.dir = viewToWorld * (imagePlane - origin);
            // ray.dir.normalize();

            // Color col = shadeRay(ray, scene, light_list, 0);
            // image.setColorAtPixel(i, j, col);
            // Averaging colors by total number of samples taken
            color[0] = color[0] / (samples*samples);
            color[1] = color[1] / (samples*samples);
            color[2] = color[2] / (samples*samples);
            color.clamp();
            image.setColorAtPixel(i, j, color);
        }
    }
}
