/***********************************************************
        

        Starter code for Assignment 3

        Implementations of functions in raytracer.h,
        and the main function which specifies the scene to be rendered.

***********************************************************/

#include "raytracer.h"
#include "util.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

const int MAX_REFLECT = 1;
const double ERR = 1e-3;

double myrand() { return rand() / ((double)RAND_MAX); }

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

// void Raytracer::computeShading(Ray3D &ray, LightList &light_list, Scene &scene) {

//     for (size_t i = 0; i < light_list.size(); ++i) {
//         LightSource *light = light_list[i];

//         // Each lightSource provides its own shading function.
//         // Implement shadows here if needed.

//         Vector3D shadow_dir = light->get_position() - ray.intersection.point;
//         shadow_dir.normalize();
//         Point3D shadow_origin = ray.intersection.point + ERR * shadow_dir;
//         // Point3D shadow_origin = ray.intersection.point;

//         Ray3D shadow_ray(shadow_origin, shadow_dir);

//         // Find intersections with ray
//         traverseScene(scene, shadow_ray);

//         if (shadow_ray.intersection.none) {
//             light->shade(ray);
//         }
//         // else {
//         //     ray.col = ray.col + ray.intersection.mat->ambient * light->get_ambient();
//         // }
//     }
// }


void Raytracer::computeShading(Ray3D &ray, LightList &light_list, Scene &scene) {
    for (size_t i = 0; i < light_list.size(); ++i) {
        LightSource *light = light_list[i];

        // Each lightSource provides its own shading function.
        // Implement shadows here if needed.
        light->shade(ray);
    }
}


Color Raytracer::shadeRay(Ray3D &ray, Scene &scene, LightList &light_list, int num_reflect) {
    Color col(0.0, 0.0, 0.0);
    traverseScene(scene, ray);

    // Don't bother shading if the ray didn't hit
    // anything.
    // if (!ray.intersection.none && (ray.intersection.mat->has_texture)) {
    //     col = textureColor(ray);
    //     return col;
    // }

    if (!ray.intersection.none) {
        computeShading(ray, light_list, scene);
        col = ray.col;
    }

    // You'll want to call shadeRay recursively (with a different ray,
    // of course) here to implement reflection/refraction effects.
    // Reflect a ray MAX_REFLECT number of times
    if (!ray.intersection.none && ray.intersection.mat->reflective && (num_reflect < MAX_REFLECT)) {
        Color reflect_color(0.0, 0.0, 0.0);

        Intersection intersect = ray.intersection;
        Vector3D normal = intersect.normal;
        Vector3D light_ray = -ray.dir;
        Vector3D reflect_dir = 2.0 * normal.dot(light_ray) * normal - light_ray;
        reflect_dir.normalize();
        Point3D reflect_origin = intersect.point + ERR * normal;

        // Our new reflected ray
        Ray3D reflect_ray(reflect_origin, reflect_dir);

        // Glossy reflection from tutorial
        if (intersect.mat->glossy) {
            int n = 10;
            double roughness = intersect.mat->roughness;
            
            // for (size_t i = 0; i < n; ++i) {
            //     Vector3D u = reflect_dir.cross(normal);
            //     u.normalize();
            //     Vector3D v = reflect_dir.cross(u);
            //     v.normalize();
            //     double theta = 2 * M_PI * ((myrand() * roughness));
            //     double phi = 2 * M_PI * ((myrand() * roughness));
            //     double x = sin(theta) * cos(phi);
            //     double y = sin(theta) * sin(phi);
            //     double z = cos(theta);

            //     reflect_dir = x * u + y * v + z * reflect_dir;
            //     reflect_dir.normalize();
            //     reflect_origin = intersect.point + ERR * reflect_dir;
            //     reflect_ray.dir = reflect_dir;
            //     reflect_ray.origin = reflect_origin;
            //     reflect_color = reflect_color + (1.0 / n) * shadeRay(reflect_ray, scene,light_list, num_reflect+1);
            // }


            Vector3D u = reflect_dir.cross(normal);
            u.normalize();
            Vector3D v = reflect_dir.cross(u);
            v.normalize();
            double theta = 2 * M_PI * ((myrand() * roughness));
            double phi = 2 * M_PI * ((myrand() * roughness));
            double x = sin(theta) * cos(phi);
            double y = sin(theta) * sin(phi);
            double z = cos(theta);

            reflect_dir = x * u + y * v + z * reflect_dir;
            reflect_dir.normalize();
            reflect_origin = intersect.point + ERR * normal;
            reflect_ray.dir = reflect_dir;
            reflect_ray.origin = reflect_origin;
            reflect_color = reflect_color + shadeRay(reflect_ray, scene, light_list, num_reflect + 1);
        } else {
            // Prevent infinite recursion
            reflect_color = shadeRay(reflect_ray, scene, light_list, num_reflect + 1);
        }

        reflect_ray.col = reflect_color;
        col = col + intersect.mat->specular * reflect_ray.col;
    }

    col.clamp();
    return col;
}

void Raytracer::renderAntiAliasDoF(Camera &camera, Scene &scene, LightList &light_list,
                                   Image &image) {
    computeTransforms(scene);
    Matrix4x4 viewToWorld;
    double factor = (double(image.height) / 2) / tan(camera.fov * M_PI / 360.0);

    viewToWorld = camera.initInvViewMatrix();
    // Construct a ray for each pixel.
    double samples = 5;
    for (int i = 0; i < image.height; i++) {
        for (int j = 0; j < image.width; j++) {
            // Implement anti aliasing with samples*samples number
            // of rays for each pixel
            Color color(0.0, 0.0, 0.0);
            for (int k = 0; k < samples; k++) {
                for (int l = 0; l < samples; l++) {
                    double rand_width = ((double)rand() / RAND_MAX) / samples + l / samples + j;
                    double rand_height = ((double)rand() / RAND_MAX) / samples + k / samples + i;

                    // Sets up ray origin and direction in view space,
                    // image plane is at z = -1.
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
                    // color = color + shadeRay(ray, scene, light_list, 0);

                    double focal_len = 10.0;
                    double apeture_size = 0.5;
                    int dof_num = 20;
                    Color dofCol = Color(0.0, 0.0, 0.0);
                    for (int count = 0; count < 3; count++) {
                        Ray3D secondary_ray =
                            computeDepthOfField(ray, origin, focal_len, apeture_size);
                        dofCol = dofCol + shadeRay(secondary_ray, scene, light_list, 0);
                    }
                    color = color + double(1.0 / dof_num) * dofCol;
                }
            }
            // Averaging colors by total number of samples taken
            //  normalize the colors and clamp them
            color[0] = color[0] / (samples * samples);
            color[1] = color[1] / (samples * samples);
            color[2] = color[2] / (samples * samples);
            color.clamp();
            image.setColorAtPixel(i, j, color);
        }
    }
}

void Raytracer::renderDoF(Camera &camera, Scene &scene, LightList &light_list, Image &image) {
    computeTransforms(scene);
    Matrix4x4 viewToWorld;
    double factor = (double(image.height) / 2) / tan(camera.fov * M_PI / 360.0);

    viewToWorld = camera.initInvViewMatrix();
    // Construct a ray for each pixel.
    for (int i = 0; i < image.height; i++) {
        for (int j = 0; j < image.width; j++) {
            // Sets up ray origin and direction in view space,
            // image plane is at z = -1.
            Color col(0.0, 0.0, 0.0);
            Point3D origin(0, 0, 0);
            Point3D imagePlane;
            imagePlane[0] = (-double(image.width) / 2 + 0.5 + j) / factor;
            imagePlane[1] = (-double(image.height) / 2 + 0.5 + i) / factor;
            imagePlane[2] = -1;

            Ray3D ray;
            // TODO: Convert ray to world space
            ray.origin = viewToWorld * origin;
            ray.dir = viewToWorld * (imagePlane - origin);
            ray.dir.normalize();
            // Color col = shadeRay(ray, scene, light_list, 0);

            // Averaging colors by total number of samples taken
            double focal_len = 3.3;
            double apeture_size = 0.2;
            int dof_num = 20;
            Color dofCol = Color(0.0, 0.0, 0.0);
            for (int count = 0; count < dof_num; count++) {
                Ray3D secondary_ray = computeDepthOfField(ray, origin, focal_len, apeture_size);
                dofCol = dofCol + shadeRay(secondary_ray, scene, light_list, 0);
            }
            col = col + double(1.0 / dof_num) * dofCol;
            col.clamp();
            image.setColorAtPixel(i, j, col);
        }
    }
}

void Raytracer::render(Camera &camera, Scene &scene, LightList &light_list, Image &image) {
    computeTransforms(scene);
    Matrix4x4 viewToWorld;
    double factor = (double(image.height) / 2) / tan(camera.fov * M_PI / 360.0);

    viewToWorld = camera.initInvViewMatrix();
    // Construct a ray for each pixel.
    for (int i = 0; i < image.height; i++) {
        for (int j = 0; j < image.width; j++) {
            // Sets up ray origin and direction in view space,
            // image plane is at z = -1.
            Point3D origin(0, 0, 0);
            Point3D imagePlane;
            imagePlane[0] = (-double(image.width) / 2 + 0.5 + j) / factor;
            imagePlane[1] = (-double(image.height) / 2 + 0.5 + i) / factor;
            imagePlane[2] = -1;

            Ray3D ray;
            // TODO: Convert ray to world space
            ray.origin = viewToWorld * origin;
            ray.dir = viewToWorld * (imagePlane - origin);
            ray.dir.normalize();
            Color col = shadeRay(ray, scene, light_list, 0);
            image.setColorAtPixel(i, j, col);
        }
    }
}

Ray3D Raytracer::computeDepthOfField(Ray3D &ray, Point3D &origin, double F, double R) {
    double normalizer = double(RAND_MAX);
    // F = focal
    // R = aperture
    // focal point
    Point3D fp = origin + F * ray.dir;

    // random point on lens
    double rand1 = R * (2 * (rand() / normalizer) - 1);
    double rand2 = R * (2 * (rand() / normalizer) - 1);
    // std::cout << rand1 << std::endl;
    // std::cout <<rand2 << std::endl;
    Point3D secondary_ray_point = origin + Vector3D(rand1, rand2, 0.0);

    // compute result
    Vector3D secondary_ray_dir = fp - secondary_ray_point;
    secondary_ray_dir.normalize();

    return Ray3D(secondary_ray_point, secondary_ray_dir);
}

Color Raytracer::textureColor(Ray3D &ray) {
    if (!ray.intersection.none && ray.intersection.mat->has_texture) {
        int width = ray.intersection.mat->width;
        int height = ray.intersection.mat->height;
        int x = ray.intersection.tex_u * width;
        int y = ray.intersection.tex_v * height;

        // Convert to rgb array index
        int i = (x * width) + y;

        double r = (ray.intersection.mat->rarray[i]) / 255.0;
        double g = (ray.intersection.mat->garray[i]) / 255.0;
        double b = (ray.intersection.mat->barray[i]) / 255.0;
        return Color(r, g, b);
    }

    return Color(0.0, 0.0, 0.0);
}
