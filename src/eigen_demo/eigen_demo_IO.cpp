#include <iostream>

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

bool parseJson(const std::string& json_whole_path)
{
  using namespace rapidjson;

  FILE* fp = fopen(json_whole_path.c_str(), "r");

  assert(fp);

  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  Document document;

  if(document.ParseStream(is).HasParseError())
  {
    std::cout << "ERROR parsing the input json file!\n";
    return false;
  }

  fclose(fp);

  assert(document.HasMember("assembly_type"));
  std::string at = document["assembly_type"].GetString();
  std::cout << "assembly type: " << at << std::endl;

  std::cout << "node list size " << document["node_list"].Size() << std::endl;

  return true;
}