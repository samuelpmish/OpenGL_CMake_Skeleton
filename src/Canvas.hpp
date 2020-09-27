#pragma once

#include <array>

#include "Shader.hpp"
#include "Application.hpp"

class Canvas : public Application {
 public:
  Canvas();

  void mouse_func(GLFWwindow* window, int button, int action, int mods);

 protected:
  virtual void loop();

 private:
  float time = 0.f;
  const int size = 100;

  // shader
  Shader vertexShader;
  Shader fragmentShader;
  ShaderProgram shaderProgram;

  int selection;
  std::array< glm::vec2, 5 > p;
  glm::mat3 W;

  // VBO/VAO/ibo
  GLuint vao, vbo, ibo;
};