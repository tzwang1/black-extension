/***********************************************************
        

        Starter code for Assignment 3

***********************************************************/

#include "raytracer.h"
#include <cstdlib>
#include "bmp_io.h"

int main(int argc, char *argv[]) {
    // Build your scene and setup your camera here, by calling
    // functions from Raytracer.  The code here sets up an example
    // scene and renders it from two different view points, DO NOT
    // change this if you're just implementing part one of the
    // assignment.
    Raytracer raytracer;
    LightList light_list;
    Scene scene;

    int width = 320;
    int height = 240;

    if (argc == 3) {
        width = atoi(argv[1]);
        height = atoi(argv[2]);
    }

    // bmp_read(char const *file_in_name, unsigned long int *width,
    //           long int *height, unsigned char **rarray, unsigned char **garray,
    //           unsigned char **barray)
    // unsigned char **rarray;
    // unsigned char **garray;
    // unsigned char **barray;

    // unsigned long floor_width = 50;
    // long floor_height = 50;

    // bool val = bmp_read("textures/horizontal.bmp", &floor_width, &floor_height, rarray, garray, barray);
    // std::cout << rarray << std::endl;

    // Define materials for shading.
    Material chrome(Color(0.25, 0.25,	0.25), Color(0.4, 0.4, 0.4),
                  Color(0.774597, 0.774597, 0.774597), 76.8);
    Material gold(Color(0.3, 0.3, 0.3), Color(0.75164, 0.60648, 0.22648),
                   Color(0.628281, 0.555802, 0.366065), 51.2);
    Material jade(Color(0, 0, 0), Color(0.54, 0.89, 0.63),
                  Color(0.316228, 0.316228, 0.316228), 12.8);

    // Defines a point light source.
    PointLight *pLight = new PointLight(Point3D(0, 0, 5), Color(0.9, 0.9, 0.9));
    light_list.push_back(pLight);

    // Add a unit sphere and square into the scene with material mat.
    SceneNode *gold_sphere = new SceneNode(new UnitSphere(), &gold);
    scene.push_back(gold_sphere);
    SceneNode *plane = new SceneNode(new UnitSquare(), &jade);
    scene.push_back(plane);

    // Add a chrome sphere
    SceneNode *chrome_sphere = new SceneNode(new UnitSphere(), &chrome);
    scene.push_back(chrome_sphere);

    // Add a cylinder
    SceneNode *gold_cylinder = new SceneNode(new UnitCylinder(), &gold);
    scene.push_back(gold_cylinder);

    // Add three walls to demonstrate reflection.
    SceneNode *wall1 = new SceneNode(new UnitSquare(), &jade);
    scene.push_back(wall1);
    wall1->mat->reflective = true;

    // SceneNode *wall2 = new SceneNode(new UnitSquare(), &jade);
    // scene.push_back(wall2);

    // SceneNode *wall3 = new SceneNode(new UnitSquare(), &jade);
    // scene.push_back(wall3);

    // Apply some transformations to the gold sphere
    double factor1[3] = {1.0, 1.0, 1.0};
    gold_sphere->translate(Vector3D(-3, 0, -5));
    gold_sphere->rotate('x', -45);
    gold_sphere->scale(Point3D(0, 0, 0), factor1);
    gold_sphere->mat->reflective = false;

    // Apply some transformations to the floor
    double factor2[3] = {10.0, 10.0, 10.0};
    plane->translate(Vector3D(0, 0, -7));
    plane->scale(Point3D(0, 0, 0), factor2);

    // Apply some transformations to the chrome sphere
    double factor3[3] = {2.0, 2.0, 2.0};
    chrome_sphere->translate(Vector3D(3, 0, -6));
    chrome_sphere->scale(Point3D(0, 0, 0), factor3);
    chrome_sphere->mat->reflective = false;

    // Apply some transformations to the gold cylinder
    gold_cylinder->translate(Vector3D(0, 0, -2));
    gold_cylinder->mat->reflective = false;
    gold_cylinder->scale(Point3D(0, 0, 0), factor1);

    // Apply some transformations to the three walls.
    wall1->translate(Vector3D(0, 3, -7));
    wall1->scale(Point3D(0, 0, 0), factor2);
    wall1->rotate('x', 45);
    // wall2->translate(Vector3D(3,0, -5));
    // wall3->translate(Vector3D(5,5,-5));

    // Render the scene, feel free to make the image smaller for
    // testing purposes.
    Camera camera1(Point3D(0, 0, 1), Vector3D(0, 0.2, -1), Vector3D(0, 1, 0),
                   80.0);
    Image image1(width, height);
    raytracer.render(camera1, scene, light_list,
                     image1);             // render 3D scene to image
    image1.flushPixelBuffer("rec_ray1.bmp"); // save rendered image to file

    // Render it from a different point of view.
    Camera camera2(Point3D(4, 2, 1), Vector3D(-4, -2, -6), Vector3D(0, 1, 0),
                   60.0);
    Image image2(width, height);
    raytracer.render(camera2, scene, light_list, image2);
    image2.flushPixelBuffer("rec_ray2.bmp");

    // Free memory
    for (size_t i = 0; i < scene.size(); ++i) {
        delete scene[i];
    }

    for (size_t i = 0; i < light_list.size(); ++i) {
        delete light_list[i];
    }

    return 0;
}