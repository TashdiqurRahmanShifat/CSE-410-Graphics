# CSE-410 Computer Graphics Course Projects

This repository contains OpenGL-based computer graphics projects developed as part of the CSE-410 Computer Graphics course. Each project demonstrates different aspects of computer graphics programming including 3D transformations, lighting, texturing, and animation.

## Prerequisites

Before running any of the projects, ensure you have the following installed:

- **GCC/MinGW compiler** with C++ support
- **OpenGL libraries**:
  - `opengl32` (usually comes with Windows)
  - `glu32` (OpenGL Utility Library)
  - `freeglut` (OpenGL Utility Toolkit)

### Installing FreeGLUT on Windows

1. Download FreeGLUT from the official website
2. Extract and place the libraries in your MinGW/lib folder
3. Place the header files in your MinGW/include/GL folder

## Compilation and Execution

All C++ files in this repository can be compiled using the following command pattern:

```bash
g++ [filename].cpp -o [output_name].exe -lglu32 -lopengl32 -lfreeglut
```

**Example:**
```bash
g++ 2005102.cpp -o rotating_fan.exe -lglu32 -lopengl32 -lfreeglut
```

Then run the executable:
```bash
./[output_name].exe
```

## Projects Overview

### 1. CSE410 OpenGL Online (B2) - Rotating Fan
**File:** `CSE410 OpenGL Online (B2)/2005102.cpp`

A 3D animated rotating fan simulation demonstrating:
- 3D object modeling with cubes and geometric shapes
- Rotation animations
- OpenGL transformations
- Camera positioning and perspective projection

**Features:**
- Interactive 3D fan with rotating blades
- Realistic 3D modeling
- Smooth animation loops

### 2. Offline 1 - Physics Simulations

#### Ball Physics Simulation
**File:** `Offline1/ball.cpp`

An advanced physics simulation featuring:
- Realistic ball physics with gravity
- Collision detection and response
- Velocity vectors and restitution
- 3D ball rotation and spin effects
- Interactive camera controls

**Controls:**
- Camera movement with keyboard
- Physics parameter adjustments
- Real-time simulation toggle

#### Analog Clock
**File:** `Offline1/clock.cpp`

A real-time analog clock implementation:
- Real-time hour, minute, and second hands
- Accurate time display
- Smooth hand movements
- Clean geometric design

### 3. Offline 2 - Graphics Pipeline Implementation
**File:** `Offline2/2005102.cpp`

A comprehensive graphics pipeline implementation including:
- 3D transformations (translation, rotation, scaling)
- Matrix operations and homogeneous coordinates
- Perspective and orthographic projections
- Clipping algorithms
- Rasterization techniques
- Bitmap image output using `bitmap_image.hpp`

**Features:**
- Complete 3D graphics pipeline from scratch
- Mathematical foundations of computer graphics
- Custom rendering without high-level OpenGL functions

### 4. Offline 3 - Ray Tracing Engine
**Files:** `Offline3/2005102.cpp`, `Offline3/2005102.h`

An advanced ray tracing renderer featuring:
- Ray-object intersection algorithms
- Multiple light sources (point lights, spot lights)
- Realistic lighting models (Phong shading)
- Texture mapping support
- Reflection and refraction effects
- Scene file parsing

**Supported Features:**
- **Geometric Objects:** Spheres, planes, triangles
- **Lighting:** Ambient, diffuse, and specular lighting
- **Textures:** Image-based texture mapping
- **Advanced Effects:** Recursive ray tracing for reflections

**Input:** Scene description file (`scene.txt`)
**Output:** High-quality rendered images

#### Bonus Implementation
**Files:** `Offline3/Bonus/2005102.cpp`, `Offline3/Bonus/2005102.h`

Enhanced ray tracing with additional advanced features.

## File Structure

```
CSE-410-Graphics/
├── README.md
├── CSE410 OpenGL Online (B2)/
│   ├── 2005102.cpp                    # Rotating fan simulation
│   ├── CSE410_Rotating_Fan.pdf        # Project specification
│   └── rotating_fan.exe               # Compiled executable
├── Offline1/
│   ├── ball.cpp                       # Ball physics simulation
│   ├── clock.cpp                      # Real-time analog clock
│   └── CSE410_Jan_25_OpenGL_Offline_Specification.pdf
├── Offline2/
│   ├── 2005102.cpp                    # Graphics pipeline implementation
│   ├── bitmap_image.hpp               # Bitmap image library
│   ├── Offline-2-Specifications.pdf  # Project specification
│   └── Test Cases.zip                 # Test cases and examples
└── Offline3/
    ├── 2005102.cpp                    # Ray tracing engine
    ├── 2005102.h                      # Header file with class definitions
    ├── bitmap_image.hpp               # Image output library
    ├── stb_image.h                    # Image loading library
    ├── scene.txt                      # Scene description file
    ├── Texture.jpg                    # Sample texture files
    ├── Texture2.jpg
    ├── Texture3.jpg
    ├── CSE410_Offline_3.pdf          # Project specification
    └── Bonus/                         # Enhanced ray tracing implementation
        ├── 2005102.cpp
        └── 2005102.h
```

## Key Learning Outcomes

Through these projects, the following computer graphics concepts are demonstrated:

1. **3D Mathematics:** Vector operations, matrix transformations, coordinate systems
2. **Rendering Pipeline:** Vertex processing, clipping, rasterization
3. **Lighting Models:** Phong lighting, multiple light sources
4. **Animation:** Frame-based animation, smooth transitions
5. **Physics Simulation:** Collision detection, realistic motion
6. **Advanced Rendering:** Ray tracing, reflections, texture mapping

## Dependencies

- **OpenGL:** Core graphics library
- **GLUT/FreeGLUT:** Window management and user input
- **STB Image:** Image loading library (for Offline 3)
- **Custom Libraries:** `bitmap_image.hpp` for image output