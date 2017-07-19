#include <iostream>

extern "C"
{
    //   ,-- 2 Leading underscores to start
    //   | ,-- then the module name
    //   | |           ,-- then _MOD_
    //   | |           |   ,-- then the subroutine name
    //   V V           V   V
    void __badlandsbmi_MOD_initialize(const char* buffer, int len);
    void __badlandsbmi_MOD_run();
    void __badlandsbmi_MOD_finalize();
}

int main() 
{

    CString fileName = "test.udm";
    const char* buffer = fileName.GetBuffer();
    int len = fileName.GetLength();
    //F_VALIDATE_XML(buffer, len);

    __badlandsbmi_MOD_initialize(buffer, len);
    __badlandsbmi_MOD_run();
    __badlandsbmi_MOD_finalize();

    return 0;
}
