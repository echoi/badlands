#include <iostream>

extern "C"
{
    //   ,-- 2 Leading underscores to start
    //   | ,-- then the module name
    //   | |           ,-- then _MOD_
    //   | |           |   ,-- then the subroutine name
    //   V V           V   V
    void __badlandsbmi_MOD_initialize();
    void __badlandsbmi_MOD_run();
    void __badlandsbmi_MOD_finalize();
}

int main() 
{

    __badlandsbmi_MOD_initialize();
    __badlandsbmi_MOD_run();
    __badlandsbmi_MOD_finalize();

    return 0;
}
