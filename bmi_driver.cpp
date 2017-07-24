#include <iostream>
#include <string.h>

extern "C"
{
    //   ,-- 2 Leading underscores to start
    //   | ,-- then the module name
    //   | |           ,-- then _MOD_
    //   | |           |   ,-- then the subroutine name
    //   V V           V   V
    void __badlands_bmi_MOD_initialize(const char* buffer, int len);
    void __badlands_bmi_MOD_run();
    void __badlands_bmi_MOD_finalize();
}

int main() 
{

    std::string fileName = "./input_detachment.xml";
    const char* buffer = fileName.c_str();
    size_t len = strlen(buffer);

    __badlands_bmi_MOD_initialize(buffer, len);
    __badlands_bmi_MOD_run();
    __badlands_bmi_MOD_finalize();

    return 0;
}
