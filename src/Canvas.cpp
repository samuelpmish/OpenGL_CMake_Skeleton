#include "Canvas.hpp"

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include <iostream>
#include <vector>
#include <random>

#include "tensor.h"

#include "glError.hpp"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

tensor<5> f(tensor<2> x) {
  return tensor<5>{2 * x(0), 2 * x(1), x(0) * x(0), 2 * x(0) * x(1), x(1) * x(1)};
}

const std::string vert_shader(R"vert(
#version 150
in vec2 position;
void main(void) {
  gl_Position = vec4(position,0.0,1.0);
}
)vert");

const std::string frag_shader(R"frag(
#version 150
out vec4 color;
uniform float radius;
uniform float width;
uniform mat3 W;
uniform vec2 p[5];
void main(void) {       

  color = vec4(1.0, 1.0, 1.0, 1.0);

  vec3 z = vec3(1.0, gl_FragCoord.xy);

  float v = abs(dot(z, W * z));
  float dvdx = length((W * z).yz);

  float r = 1.0e10;
  for (int i = 0; i < 5; i++) {
    r = min(r, length(gl_FragCoord.xy - p[i]));
  }

  float q = min(smoothstep(radius, radius + 2, r), smoothstep(width, width + 2, v / dvdx));

  color = vec4(0.0, 0.0, 1.0 - q, 1.0);

}
)frag");

void mouse_callback(GLFWwindow* window, int button, int action, int mods) {
  auto canvas = (Canvas *)glfwGetWindowUserPointer(window);
  canvas->mouse_func(window, button, action, mods);
}

void Canvas::mouse_func(GLFWwindow* window, int button, int action, int mods) {
  std::cout << "mouse func" << std::endl;
  if (button == GLFW_MOUSE_BUTTON_LEFT) {
    if(GLFW_PRESS == action) {
      double x, y;
      glfwGetCursorPos(window, &x, &y);
      glm::vec2 m{float(x), float(height - y)};

      float distance = 100.0f;
      for (int i = 0; i < 5; i++) {
        if (glm::distance(p[i], m) < distance) {
          selection = i;
          distance = glm::distance(p[i], m);
        }
      }
    } else {
      selection = -1;
    }
  }
}

Canvas::Canvas() : Application(),
      vertexShader(Shader::fromString(vert_shader, GL_VERTEX_SHADER)),
      fragmentShader(Shader::fromString(frag_shader, GL_FRAGMENT_SHADER)),
      shaderProgram({vertexShader, fragmentShader}) {
  glCheckError(__FILE__, __LINE__);

  selection = -1;

  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(0.0, 1.0);

  for (int i = 0; i < 5; i++) {
    p[i][0] = distribution(generator) * 1000.0f;
    p[i][1] = distribution(generator) * 1000.0f;
  }

  // creation of the mesh ------------------------------------------------------
  std::vector<glm::vec2> vertices;

  vertices.push_back({-1.0, -1.0});
  vertices.push_back({ 1.0, -1.0});
  vertices.push_back({ 1.0,  1.0});

  vertices.push_back({ 1.0,  1.0});
  vertices.push_back({-1.0,  1.0});
  vertices.push_back({-1.0, -1.0});

  // creation of the vertex array buffer----------------------------------------

  // vbo
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec2),
               vertices.data(), GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  // vao
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  // bind vbo
  glBindBuffer(GL_ARRAY_BUFFER, vbo);

  // map vbo to shader attributes
  shaderProgram.setAttribute("position", 2, sizeof(glm::vec2), 0);

  // vao end
  glBindVertexArray(0);

  glfwSetWindowUserPointer(window, (void *)this);
  glfwSetMouseButtonCallback(window, mouse_callback);

}

void Canvas::loop() {

  static float width = 3.0;
  static float radius = 4.0;

  // exit on window close button pressed
  if (glfwWindowShouldClose(getWindow()))
    exit();

  float t = getTime();

  // clear
  glClear(GL_COLOR_BUFFER_BIT);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // feed inputs to dear imgui, start new frame
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  shaderProgram.use();

  // send uniforms
  shaderProgram.setUniform("p[0]", p[0]);
  shaderProgram.setUniform("p[1]", p[1]);
  shaderProgram.setUniform("p[2]", p[2]);
  shaderProgram.setUniform("p[3]", p[3]);
  shaderProgram.setUniform("p[4]", p[4]);

  if (selection != -1) {
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    p[selection][0] = float(x);
    p[selection][1] = float(height - y);

    tensor<5,5> A{};
    tensor<5> b{};
    for (int i = 0; i < 5; i++) {
      auto f0 = f({p[i][0], p[i][1]});
      A += outer(f0, f0);
      b += -f0;
    }
    auto z = linear_solve(A, b);

    W[0][0] = 1.0f;
    W[1][0] = W[0][1] = float(z(0));
    W[2][0] = W[0][2] = float(z(1));
    W[1][1] = float(z(2));
    W[2][1] = W[1][2] = float(z(3));
    W[2][2] = float(z(4));

    std::cout << (W[2][2] * W[1][1] - W[1][2] * W[2][1]) << std::endl;
  }

  shaderProgram.setUniform("W", W);
  shaderProgram.setUniform("radius", radius);
  shaderProgram.setUniform("width", width);

  glCheckError(__FILE__, __LINE__);

  glBindVertexArray(vao);

  glBindBuffer(GL_ARRAY_BUFFER, vbo);

  glCheckError(__FILE__, __LINE__);
  glDrawArrays(GL_TRIANGLES, 0, 6);

  glBindVertexArray(0);

  shaderProgram.unuse();

  // render your GUI
  ImGui::Begin("Demo window");

  ImGui::DragFloat("radius", &radius, 0.1f, 1.0, 15.0);
  ImGui::DragFloat("width", &width, 0.1f, 1.0, 15.0);

  ImGui::End();

  // Render dear imgui into screen
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

}
