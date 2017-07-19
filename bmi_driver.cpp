extern "C"
{
    //          ,-- 2 Leading underscores to start
    //          | ,-- then the module name
    //          | |            ,-- then _MOD_
    //          | |            |    ,-- then the subroutine name
    //          V V            V    V
    extern void __BADLANDS_BMI_MOD_BADLANDS_Initialize();
    extern void __BADLANDS_BMI_MOD_BADLANDS_Run();
    extern void __BADLANDS_BMI_MOD_BADLANDS_Finalize();
}

int main() 
{

    __BADLANDS_BMI_MOD_BADLANDS_Initialize();
    __BADLANDS_BMI_MOD_BADLANDS_Run();
    __BADLANDS_BMI_MOD_BADLANDS_Finalize();

}
