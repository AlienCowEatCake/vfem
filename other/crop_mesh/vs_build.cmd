@call "%VS71COMNTOOLS%\vsvars32.bat"
@cd %~dp0
cl /c /nologo /Ox /EHsc /W3 /MT /DNDEBUG crop_mesh.cpp /Focrop_mesh.cpp.obj
@if not exist "crop_mesh.cpp.obj" goto err
link /NOLOGO /SUBSYSTEM:CONSOLE crop_mesh.cpp.obj /OUT:crop_mesh.exe
@if not exist "crop_mesh.exe" goto err
@goto end
:err
@echo Error!
:end
@pause
