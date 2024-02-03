# NoopMath

NoopMath is a bare bones linear algebra library with emphasis on integration with OpenGL.

## Specifics
NoopMath stores matrices column-major. So, for a Mat4, the first 4 elements make up the first
column of the matrix. For a Mat4, elements at index 0, 4, 8, and 12 make up the first row of the
matrix. Vectors are treated as column vectors. This decision is strictly to maintain the standard
conventions of mathematics and OpenGL itself.

NoopMath is agnostic to left-hand/right-hand. It's just a linear algebra number cruncher.
Multiply two matrices and it will do matrix multiplication. As simple as that.

## Use

### Create a Camera Matrix

    // camera.origin is a vec3
    mat4 translation{
                     1.0f,              0.0f,              0.0f, 0.0f,
                     0.0f,              1.0f,              0.0f, 0.0f,
                     0.0f,              0.0f,              1.0f, 0.0f,
        -camera.origin[0], -camera.origin[1], -camera.origin[2], 1.0f
    };

    // camera.right, camera.up, and camera.forward are all vec3
    mat4 measure{
        camera.right[0], camera.up[0], -camera.forward[0], 0.0f,
        camera.right[1], camera.up[1], -camera.forward[1], 0.0f,
        camera.right[2], camera.up[2], -camera.forward[2], 0.0f,
                   0.0f,         0.0f,               0.0f, 1.0f
    };

    // Matrices should be read from right to left
    // translate object, so that the camera lies at the origin
    // measure the position of the object in the direction the camera is facing
    mat4 cameraMat = measure * translation;

### Uniform Helper Functions

    inline void setUniform(GLuint shaderId, const std::string& name, const mat4* mat) {
        glUniformMatrix4fv(glGetUniformLocation(shaderId, name.c_str()), 1, GL_FALSE, mat->values); 
    }

    inline void setUniform(GLuint shaderId, const std::string& name, const vec2& vector2) {
        glUniform2fv(glGetUniformLocation(shaderId, name.c_str()), 1, vector2.values);
    }

    inline void setUniform(GLuint shaderId, const std::string& name, const vec3& vector3) {
        glUniform3fv(glGetUniformLocation(shaderId, name.c_str()), 1, vector3.values);
    }

    inline void setUniform(GLuint shaderId, const std::string& name, const vec4& vector4) {
        glUniform4fv(glGetUniformLocation(shaderId, name.c_str()), 1, vector4.values);
    }
