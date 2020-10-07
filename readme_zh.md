
## 网格化结果保存为stl，供snappyhexmesh使用

### 开发流程

1. 从Fltk/graphicWindow.cpp源文件入手，找到`_save_stl`函数，然后定位到`stlFileDialog`函数

2. `stlFileDialog`函数位于Fltk/fileDialogs.cpp原文件中，在此函数最后面有一个`CreateOutputFile`函数，这是真正实现写stl的函数

3. `CreateOutputFile`函数位于Common/CreateFile.cpp源文件中，这里有switch...case，找到`case FORMAT_STL`，里面有一个`writeSTL`函数，这就是关键

4. `writeSTL`函数位于Geo/GModelIO_STL.cpp文件中，在此函数中就可以将网格以physical分组进行分别写入Surface和Volume的stl文件，在stl文件里面用physical名称命名不同的solid。相比修改代码前，新增了两个stl文件，一个用xxx_Surfaces.stl命名，另一个用xxx_Volumes.stl命名（如果有volume的话）


