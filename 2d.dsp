# Microsoft Developer Studio Project File - Name="2d" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=2d - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "2d.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "2d.mak" CFG="2d - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "2d - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "2d - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "2d - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "2d - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /fpscomp:general /fpscomp:logicals /fpscomp:ioformat /nologo /traceback /warn:argument_checking /warn:nofileopt
# SUBTRACT F90 /fpscomp:ldio_spacing
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "2d - Win32 Release"
# Name "2d - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\arrayDeinitialisation.f90
DEP_F90_ARRAY=\
	".\Debug\commonArrays.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\arrayInitialisation.f90
DEP_F90_ARRAYI=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\calculation2D.f90
DEP_F90_CALCU=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\checkPrint.f90
DEP_F90_CHECK=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\commonArrays.f90
# End Source File
# Begin Source File

SOURCE=.\commonVariables.f90
# End Source File
# Begin Source File

SOURCE=.\conditionsBound.f90
DEP_F90_CONDI=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\conditionsInlet.f90
DEP_F90_CONDIT=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\conditionsOutlet.f90
DEP_F90_CONDITI=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\crash.f90
DEP_F90_CRASH=\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\derivatives.f90
DEP_F90_DERIV=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\eigenvalues.f90
DEP_F90_EIGEN=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\fl_X.f90
DEP_F90_FL_X_=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\fl_Y.f90
DEP_F90_FL_Y_=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\fluxCalc.f90
DEP_F90_FLUXC=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\geom.f90
DEP_F90_GEOM_=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\initial.f90
DEP_F90_INITI=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\isnas.f90
# End Source File
# Begin Source File

SOURCE=.\main.f90
DEP_F90_MAIN_=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\minmod.f90
# End Source File
# Begin Source File

SOURCE=.\multAeta.f90
DEP_F90_MULTA=\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\multAksi.f90
DEP_F90_MULTAK=\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\multL.f90
DEP_F90_MULTL=\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\multLinv.f90
DEP_F90_MULTLI=\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\MultT.f90
DEP_F90_MULTT=\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PrintFromF.f90
DEP_F90_PRINT=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\readFiles.f90
DEP_F90_READF=\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\residual.f90
DEP_F90_RESID=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\RR1d.f90
DEP_F90_RR1D_=\
	".\Debug\commonVariables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\timeStep.f90
DEP_F90_TIMES=\
	".\Debug\commonArrays.mod"\
	".\Debug\commonVariables.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
