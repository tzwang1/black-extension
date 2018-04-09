/***********************************************************
        
    Starter code for Assignment 3

    Implements scene_object.h

    Equations for primitives taken from
    https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html#SECTION00023200000000000000

***********************************************************/

#include "scene_object.h"
#include <cmath>

// Checks whether values for a, b, and c solve a quadratic equation
bool solveQuadratic(double &a, double &b, double &c, double &t) {
    double discr = b * b - 4 * a * c;
    double x0, x1;

    if (discr < 0) {
        return false;
    } else if (discr == 0) {
        x0 = -b / (2 * a);
        x1 = -b / (2 * a);
    } else {
        x0 = (-b + sqrt(discr)) / (2 * a);
        x1 = (-b - sqrt(discr)) / (2 * a);
    }
    // solve for parameter t by taking the minimum between the two solns
    // return false if the t is negative
    if (x0 < 0 && x1 < 0) {
        return false;
    } else if (x0 < 0) {
        t = x1;
    } else if (x1 < 0) {
        t = x0;
    } else {
        t = fmin(x0, x1);
    }

    return true;
}

// Sets values for an intersecting point
void set_intersect_values(Ray3D &ray, double t, Vector3D normal, Point3D point,
                          const Matrix4x4 &worldToModel, const Matrix4x4 &modelToWorld) {
    ray.intersection.t_value = t;
    ray.intersection.point = modelToWorld * point;
    normal = transNorm(worldToModel, normal);
    normal.normalize();
    ray.intersection.normal = normal;
    ray.intersection.none = false;
}

bool UnitSquare::intersect(Ray3D &ray, const Matrix4x4 &worldToModel,
                           const Matrix4x4 &modelToWorld) {
    // TODO: implement intersection code for UnitSquare, which is
    // defined on the xy-plane, with vertices (0.5, 0.5, 0),
    // (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
    // (0, 0, 1).
    //
    // Your goal here is to fill ray.intersection with correct values
    // should an intersection occur.  This includes intersection.point,
    // intersection.normal, intersection.none, intersection.t_value.
    //
    // HINT: Remember to first transform the ray into object space
    // to simplify the intersection test.

    // Transform the ray (origin, direction) to object space
    Point3D origin = worldToModel * ray.origin;
    Vector3D direction = worldToModel * ray.dir;

    // Since the normal of the unit square is (0, 0, 1),
    // we only need to consider the z coordinate.
    double t = -origin[2] / direction[2];

    // Invalid intersection
    if (t <= 0) {
        return false;
    }

    // Calculate the point where the ray intersects the xy-plane.
    Point3D point = origin + t * direction;
    Vector3D normal(0.0, 0.0, 1.0);

    // Check if the point of intersection is within the unit square.
    if ((-0.5 <= point[0]) && (point[0] <= 0.5) && (-0.5 <= point[1]) && (point[1] <= 0.5)) {
        if (ray.intersection.none || t < ray.intersection.t_value) {
            set_intersect_values(ray, t, normal, point, worldToModel, modelToWorld);
            return true;
        }
    }

    return false;
}

bool UnitSphere::intersect(Ray3D &ray, const Matrix4x4 &worldToModel,
                           const Matrix4x4 &modelToWorld) {
    // TODO: implement intersection code for UnitSphere, which is centred
    // on the origin.
    //
    // Your goal here is to fill ray.intersection with correct values
    // should an intersection occur.  This includes intersection.point,
    // intersection.normal, intersection.none, intersection.t_value.
    //
    // HINT: Remember to first transform the ray into object space
    // to simplify the intersection test.

    // t is the parameter for the point that intersects the sphere
    double t;
    // Transform the ray (origin, direction) to object space
    Point3D origin = worldToModel * ray.origin;
    Vector3D direction = worldToModel * ray.dir;

    Vector3D center = origin - Point3D(0, 0, 0);
    double a = direction.dot(direction);
    double b = 2 * direction.dot(center);
    double c = center.dot(center) - 1;

    if (!solveQuadratic(a, b, c, t)) {
        return false;
    }

    // Intersection point
    Point3D point = origin + t * direction;

    Vector3D normal(point[0], point[1], point[2]);
    normal.normalize();

    if (ray.intersection.none || t < ray.intersection.t_value) {
        set_intersect_values(ray, t, normal, point, worldToModel, modelToWorld);
        return true;
    }

    return false;
}

bool intersects_top_or_bot(double t, Point3D origin, Vector3D direction,
                           Intersection intersection) {
    Point3D point = origin + t * direction;
    if (point[0] * point[0] + point[2] * point[2] < 1 && t > 0 && t < intersection.t_value) {
        return true;
    }
    return false;
}

bool UnitCylinder::intersect(Ray3D &ray, const Matrix4x4 &worldToModel,
                             const Matrix4x4 &modelToWorld) {
    // convert ray to model space
    Ray3D rayModelSpace = Ray3D(worldToModel * ray.origin, worldToModel * ray.dir);

    // center of unit cyclinder
    Point3D center(0, 0, 0);

    double t_value = 0.0;
    double temp1, temp2;

    // the radius of the disk
    double radius = 1.0;

    Vector3D normal;

    Point3D intersectionPoint;

    if (rayModelSpace.dir[2] != 0) {
        //  intersection with z = 0.5
        temp1 = (-0.5 - rayModelSpace.origin[2]) / rayModelSpace.dir[2];
        //  intersection with z = 0.5
        temp2 = (0.5 - rayModelSpace.origin[2]) / rayModelSpace.dir[2];

        // intersect z = 0.5 first
        t_value = temp2;
        Point3D normal_temp(0, 0, 1);

        if (temp1 < temp2) {
            t_value = temp1;
            // Construct the normal for the closer disk, which points on the
            // negative z axis
            Point3D normal_temp(0, 0, -1);
        }
        normal = normal_temp - center;
        normal.normalize();

        if (t_value <= 0.0001) {

            return false;
        }

        intersectionPoint = rayModelSpace.origin + t_value * rayModelSpace.dir;

        // intersects with top/ bot disk
        if (intersectionPoint[0] * intersectionPoint[0] +
                intersectionPoint[1] * intersectionPoint[1] <=
            (radius * radius)) {

            if (ray.intersection.none || t_value < ray.intersection.t_value) {
                ray.intersection.point = intersectionPoint;

                ray.intersection.normal = normal;
                ray.intersection.t_value = t_value;
                ray.intersection.none = false;

                return true;
            }

            return false;
        }

        // intersection with body
        double a = rayModelSpace.dir[0] * rayModelSpace.dir[0] +
                   rayModelSpace.dir[1] * rayModelSpace.dir[1];
        double b = (rayModelSpace.origin[0] * rayModelSpace.dir[0] +
                    rayModelSpace.origin[1] * rayModelSpace.dir[1]);
        double c = rayModelSpace.origin[0] * rayModelSpace.origin[0] +
                   rayModelSpace.origin[1] * rayModelSpace.origin[1] - radius * radius;

        double discriminant = b * b - a * c;

        if (discriminant < 0) {
            // no real root
            return false;
        }

        double root1 = (-b + sqrt(discriminant)) / a;
        double root2 = (-b - sqrt(discriminant)) / a;

        if (root1 > 0 || root2 > 0) {

            if (root1 > 0 && root2 > 0) {

                if (root1 <= root2) {
                    t_value = root1;
                }

                else {
                    t_value = root2;
                }
            }

            else if (root1 > 0) {
                t_value = root1;
            }

            else {
                t_value = root2;
            }
            intersectionPoint = rayModelSpace.origin + t_value * rayModelSpace.dir;
            normal[0] = intersectionPoint[0];
            normal[1] = intersectionPoint[1];
            normal[2] = 0;
            normal.normalize();

            if (intersectionPoint[2] < 0.5 && intersectionPoint[2] > -0.5) {

                if (ray.intersection.none || t_value < ray.intersection.t_value) {
                    ray.intersection.point = modelToWorld * intersectionPoint;
                    Point3D normalTemp;
                    normalTemp[0] = intersectionPoint[0];
                    normalTemp[1] = intersectionPoint[1];
                    normalTemp[2] = 0;
                    ray.intersection.normal = modelToWorld * (normalTemp - center);
                    ray.intersection.t_value = t_value;
                    ray.intersection.none = false;
                    return true;
                }

                return false;
            }

            return false;
        }

        return false;
    }

    // rayModelSpace.origin[2] == 0

    double a = pow(rayModelSpace.dir[0], 2) + pow(rayModelSpace.dir[1], 2);
    double b = 2 * (rayModelSpace.dir[0] * rayModelSpace.origin[0] +
                    rayModelSpace.dir[1] * rayModelSpace.origin[1]);
    double c = pow(rayModelSpace.origin[0], 2) + pow(rayModelSpace.origin[1], 2) - 1;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {

        return false;
    }

    double root1 = (-b + sqrt(discriminant)) / (2 * a);
    double root2 = (-b - sqrt(discriminant)) / (2 * a);

    if (root1 > 0 || root2 > 0) {

        if (root1 > 0 && root2 > 0) {

            if (root1 <= root2) {
                t_value = root1;
            }

            else {
                t_value = root2;
            }
        }

        else if (root1 > 0) {
            t_value = root1;
        }

        else {
            t_value = root2;
        }

        intersectionPoint = rayModelSpace.origin + t_value * rayModelSpace.dir;
        normal[0] = intersectionPoint[0];
        normal[1] = intersectionPoint[1];
        normal[2] = 0;
        normal.normalize();

        if (intersectionPoint[2] < 0.5 && intersectionPoint[2] > -0.5) {

            if (ray.intersection.none || t_value < ray.intersection.t_value) {
                ray.intersection.point = modelToWorld * intersectionPoint;
                Point3D normalTemp;
                normalTemp[0] = intersectionPoint[0];
                normalTemp[1] = intersectionPoint[1];
                normalTemp[2] = 0;
                ray.intersection.normal = modelToWorld * (normalTemp - center);
                ray.intersection.t_value = t_value;
                ray.intersection.none = false;
                return true;
            }

            return false;
        }

        return false;
    }

    return false;
}

void SceneNode::rotate(char axis, double angle) {
    Matrix4x4 rotation;
    double toRadian = 2 * M_PI / 360.0;
    int i;

    for (i = 0; i < 2; i++) {
        switch (axis) {
        case 'x':
            rotation[0][0] = 1;
            rotation[1][1] = cos(angle * toRadian);
            rotation[1][2] = -sin(angle * toRadian);
            rotation[2][1] = sin(angle * toRadian);
            rotation[2][2] = cos(angle * toRadian);
            rotation[3][3] = 1;
            break;
        case 'y':
            rotation[0][0] = cos(angle * toRadian);
            rotation[0][2] = sin(angle * toRadian);
            rotation[1][1] = 1;
            rotation[2][0] = -sin(angle * toRadian);
            rotation[2][2] = cos(angle * toRadian);
            rotation[3][3] = 1;
            break;
        case 'z':
            rotation[0][0] = cos(angle * toRadian);
            rotation[0][1] = -sin(angle * toRadian);
            rotation[1][0] = sin(angle * toRadian);
            rotation[1][1] = cos(angle * toRadian);
            rotation[2][2] = 1;
            rotation[3][3] = 1;
            break;
        }
        if (i == 0) {
            this->trans = this->trans * rotation;
            angle = -angle;
        } else {
            this->invtrans = rotation * this->invtrans;
        }
    }
}

void SceneNode::translate(Vector3D trans) {
    Matrix4x4 translation;

    translation[0][3] = trans[0];
    translation[1][3] = trans[1];
    translation[2][3] = trans[2];
    this->trans = this->trans * translation;
    translation[0][3] = -trans[0];
    translation[1][3] = -trans[1];
    translation[2][3] = -trans[2];
    this->invtrans = translation * this->invtrans;
}

void SceneNode::scale(Point3D origin, double factor[3]) {
    Matrix4x4 scale;

    scale[0][0] = factor[0];
    scale[0][3] = origin[0] - factor[0] * origin[0];
    scale[1][1] = factor[1];
    scale[1][3] = origin[1] - factor[1] * origin[1];
    scale[2][2] = factor[2];
    scale[2][3] = origin[2] - factor[2] * origin[2];
    this->trans = this->trans * scale;
    scale[0][0] = 1 / factor[0];
    scale[0][3] = origin[0] - 1 / factor[0] * origin[0];
    scale[1][1] = 1 / factor[1];
    scale[1][3] = origin[1] - 1 / factor[1] * origin[1];
    scale[2][2] = 1 / factor[2];
    scale[2][3] = origin[2] - 1 / factor[2] * origin[2];
    this->invtrans = scale * this->invtrans;
}
