#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"



int main(void)
{

    const int screenWidth = 1024;
    const int screenHeight = 720;

    InitWindow(screenWidth, screenHeight, "SMEBotic v0.1 by Luis Santos aka DJOKER");
    SetTargetFPS(60);
 
  

    Camera camera = {0};
    camera.position = (Vector3){10.0f, 10.0f, 10.0f};
    camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    bool setCamera = false;

    while (!WindowShouldClose())
    {
        if (setCamera)
            UpdateCamera(&camera, CAMERA_FREE);

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
        {
            if (setCamera)
                setCamera = false;
            else
                setCamera = true;
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
        {
            setCamera = false;
        }
        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode3D(camera);

        DrawGrid(10, 1.0f);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){1.0f, 0.5f, 0.0f}, RED);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){0.0f, 0.5f, 1.0f}, GREEN);
        DrawLine3D((Vector3){0.0f, 0.5f, 0.0f}, (Vector3){0.0f, 1.5f, 0.0f}, BLUE);
       
        EndMode3D();

        EndDrawing();
    }

    

    CloseWindow();

    return 0;
}