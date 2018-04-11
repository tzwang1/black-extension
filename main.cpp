/***********************************************************
        

        Starter code for Assignment 3

***********************************************************/

#include "raytracer.h"
#include <cstdlib>
#include "bmp_io.h"

void renderTextureScene(int width, int height) {
    std::cout << "Rendering textured scene" <<std::endl;
    Scene scene;
    Raytracer raytracer;
    LightList light_list;

    Material floor_texture;
    if(!bmp_read("textures/tube.bmp", &floor_texture.width, &floor_texture.height, &floor_texture.rarray, &floor_texture.garray, &floor_texture.barray)){
        std::cout << "Floor texture was sucessfully read!\n" << std::endl;
    } else {
        std::cout << "Something went wrong while reading the texture" << std::endl;
    }
    floor_texture.has_texture = true;

    Material foil;
    if(!bmp_read("textures/Foil1.bmp", &foil.width, &foil.height, &foil.rarray, &foil.garray, &foil.barray)){
        std::cout << "Foil texture was successfully read!\n" << std::endl;
    } else {
        std::cout << "Something went wrong while reading the texture" << std::endl;
    }
    foil.has_texture = true;

    // Defines a point light source.
    PointLight *pLight = new PointLight(Point3D(0, 0, 5), Color(0.9, 0.9, 0.9));
    light_list.push_back(pLight);

    // Add a foil sphere
    SceneNode *foil_sphere = new SceneNode(new UnitSphere(), &foil);
    scene.push_back(foil_sphere);

    // Add a cube
    SceneNode *copper_cube = new SceneNode(new UnitCube(), &foil);
    scene.push_back(copper_cube);

    // Add a plane
    SceneNode *plane = new SceneNode(new UnitSquare(), &floor_texture);
    scene.push_back(plane);

    double factor1[3] = {1.0, 1.0, 1.0};
    double factor2[3] = {20.0, 20.0, 10.0};
    double factor3[3] = {2.0, 2.0, 2.0};

    plane->translate(Vector3D(0, 0, -7));
    plane->scale(Point3D(0, 0, 0), factor2);
    plane->mat->reflective = true;

    // Apply some transformations to the foil sphere
    foil_sphere->translate(Vector3D(3, 2, -6));
    foil_sphere->scale(Point3D(0, 0, 0), factor3);
    foil_sphere->mat->reflective = false;
    foil_sphere->mat->transparency = 0.5;
    foil_sphere->mat->refractive = true;

    // Apply some transformations to the cube
    copper_cube->translate(Vector3D(-3, -2, -4));
    copper_cube->scale(Point3D(0, 0, 0), factor3);

    // Render the scene, feel free to make the image smaller for
    // testing purposes.
    Camera camera1(Point3D(0, 0, 3), Vector3D(0, 0.2, -1), Vector3D(0, 1, 0),
                   80.0);
    Image image1(width, height);
    raytracer.render(camera1, scene, light_list,
                     image1);             // render 3D scene to image
    image1.flushPixelBuffer("image1.bmp"); // save rendered image to file

    // Render it from a different point of view.
    Camera camera2(Point3D(4, 2, 1), Vector3D(-4, -2, -7), Vector3D(0, 1, 0),
                   60.0);
    Image image2(width, height);
    raytracer.render(camera2, scene, light_list, image2);
    image2.flushPixelBuffer("image2.bmp");

    Camera camera3(Point3D(0,-5,8), Vector3D(0, 1, -1), Vector3D(0, 1, 0),
                    100.0);
    Image image3(width, height);
    raytracer.renderDoF(camera3, scene, light_list, image3);
    image3.flushPixelBuffer("image3.bmp");

    Camera camera4(Point3D(2,3,1), Vector3D(-1, -2, -1), Vector3D(0, 1, 0),
                    70.0);
    Image image4(width, height);
    raytracer.renderAntiAliasDoF(camera4, scene, light_list, image4);
    image4.flushPixelBuffer("image4.bmp");

    // Free memory
    for (size_t i = 0; i < scene.size(); ++i) {
        delete scene[i];
    }

    for (size_t i = 0; i < light_list.size(); ++i) {
        delete light_list[i];
    }
}

void renderNonTexturedScene(int width, int height){
    std::cout << "Rendering non textured scene" <<std::endl;
    Raytracer raytracer;
    LightList light_list;
    Scene scene;

    // Define materials for shading.
    Material chrome(Color(0.25, 0.25, 0.25), Color(0.4, 0.4, 0.4),
                  Color(0.774597, 0.774597, 0.774597), 76.8);
    Material gold(Color(0.3, 0.3, 0.3), Color(0.75164, 0.60648, 0.22648),
                   Color(0.628281, 0.555802, 0.366065), 51.2);
    Material jade(Color(0, 0, 0), Color(0.54, 0.89, 0.63),
                  Color(0.316228, 0.316228, 0.316228), 12.8);
    
    Material copper(Color(0.19125, 0.0735, 0.0225), Color(0.7038, 0.27048, 0.0828),
                  Color(0.256777, 0.137622, 0.086014), 11.264);

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

    // // Add a cylinder
    SceneNode *copper_cylinder = new SceneNode(new UnitCylinder(), &copper);
    scene.push_back(copper_cylinder);

    SceneNode *wall1 = new SceneNode(new UnitSquare(), &jade);
    scene.push_back(wall1);

    // Apply some transformations to the gold sphere
    double factor1[3] = {1.0, 1.0, 1.0};
    gold_sphere->translate(Vector3D(-4, 2, -5));
    gold_sphere->rotate('x', -45);
    gold_sphere->scale(Point3D(0, 0, 0), factor1);
    gold_sphere->mat->reflective = false;
    gold_sphere->mat->glossy = true;

    // Apply some transformations to the floor
    double factor2[3] = {20.0, 20.0, 10.0};
    plane->translate(Vector3D(0, 0, -7));
    plane->scale(Point3D(0, 0, 0), factor2);
    plane->mat->reflective = true;

    // Apply some transformations to the chrome sphere
    double factor3[3] = {2.0, 2.0, 2.0};
    chrome_sphere->translate(Vector3D(1, 0, -6));
    chrome_sphere->scale(Point3D(0, 0, 0), factor1);
    chrome_sphere->mat->reflective = true;
    chrome_sphere->mat->reflective_index = 1.0;
    chrome_sphere->mat->refractive = true;
    chrome_sphere->mat->refractive_index = 1.3;
    chrome_sphere->mat->glossy = true;
    chrome_sphere->mat->transparency = 0.5;
    
    // Apply some transformations to the gold cylinder
    copper_cylinder->translate(Vector3D(-2, -2, -6));
    copper_cylinder->mat->reflective = false;
    copper_cylinder->mat->glossy = true;
    copper_cylinder->scale(Point3D(0, 0, 0), factor1);

    // Apply some transformations to the wall
    double factor4[3] = {5.0, 5.0, 5.0};
    wall1->rotate('x', 45);
    wall1->translate(Vector3D(0, 5, -6));
    wall1->scale(Point3D(0,0,0), factor2);
    wall1->mat->reflective = true;
    wall1->mat->glossy = false;

    // Render the scene, feel free to make the image smaller for
    // testing purposes.
    Camera camera1(Point3D(0, 0, 3), Vector3D(0, 0.2, -1), Vector3D(0, 1, 0),
                   80.0);
    Image image1(width, height);
    raytracer.render(camera1, scene, light_list,
                     image1);             // render 3D scene to image
    image1.flushPixelBuffer("image1.bmp"); // save rendered image to file

    // Render it from a different point of view.
    Camera camera2(Point3D(4, 2, 1), Vector3D(-4, -2, -7), Vector3D(0, 1, 0),
                   60.0);
    Image image2(width, height);
    raytracer.render(camera2, scene, light_list, image2);
    image2.flushPixelBuffer("image2.bmp");

    // Camera camera3(Point3D(7,0,10), Vector3D(-1, 2, -5), Vector3D(0, 1, 0),
    //                 90.0);
    // Image image3(width, height);
    // raytracer.renderDoF(camera3, scene, light_list, image3);
    // image3.flushPixelBuffer("image3.bmp");

    // Camera camera4(Point3D(2,3,2), Vector3D(-1, -2, -1), Vector3D(0, 1, 0),
    //                 70.0);
    // Image image4(width, height);
    // raytracer.renderAntiAliasDoF(camera3, scene, light_list, image4);
    // image4.flushPixelBuffer("image4.bmp");

    // Free memory
    for (size_t i = 0; i < scene.size(); ++i) {
        delete scene[i];
    }

    for (size_t i = 0; i < light_list.size(); ++i) {
        delete light_list[i];
    }
}

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
    renderNonTexturedScene(width, height);
    // renderTextureScene(width, height);
    return 0;
}