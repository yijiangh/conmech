//#include "eigen_demo/eigen_demo_IO.cpp"

#include <string>

class Pet
{
 public:
  Pet(const std::string &name) : name(name) { }

  void setName(const std::string &name_) { name = name_; }
  const std::string &getName() const { return name; }

 private:
  std::string name;

 public:
  std::string file_path;
};