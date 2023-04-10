#include <iostream>
#include <window.hpp>
#include <qapplication.h>

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    appwindow wnd;
    wnd.setMinimumHeight(900);
    wnd.show();
    return app.exec();
}