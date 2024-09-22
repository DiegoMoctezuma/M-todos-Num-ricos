//Progama uno de la materia Métodos Numéricos
//Autores:
//  - Moctezuma Ramirez Diego Rafael
//  - Segura Loera Carlos Emiliano
//  - Ruiz Garcia Emiliano
//  - Moreno Vigueras Arturo Tadeo

#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include<math.h>

#include <ftxui/component/captured_mouse.hpp>
#include <ftxui/component/component.hpp>
#include <ftxui/component/component_options.hpp>
#include <ftxui/component/component_base.hpp>
#include <ftxui/component/screen_interactive.hpp>

#include <ftxui/screen/color.hpp>
#include <ftxui/dom/elements.hpp>
#include <ftxui/dom/flexbox_config.hpp>

using uint = unsigned int;
using namespace std;
using namespace ftxui;

struct Tabla {
    //Atributos
    vector<vector<string>> data;
    int rows,columns;
    string raiz, iteracion; 

    // Constructor para inicializar la tabla
    Tabla(int rows, int cols) {
        this->rows = rows;
        this->columns = cols;

        // Inicializar la tabla con datos de ejemplo
        for (int i = 0; i < rows; ++i) {
            vector<string> row;
            for (int j = 0; j < cols; ++j) {
                row.push_back("-");
            }
            data.push_back(row);
        }
    }

    //Llenado de la tabla dependiendo de la funcion y el metodo
    void llenarTablaBiseccion(unsigned int funcion, double a, double b, double tol, int iter){
        //Limpiar la tabla
        data.clear();
        vector<string> headers = {"i", "a", "b", "f(a)", "f(b)", "p", "f(p)", "Er"};
        data.push_back(headers);

        //Funciones auxiliares
        long double fa = 0, fb = 0, p = 0, fp = 0, Er = 0, aux = 0;
        string criterio;

        //Llenar la tabla
        for (int i = 0; i < iter + 1; ++i) {
            //Verifica la funcion
            switch(funcion){
                case 0:
                    fa = powl(a,2) * cosl(a) - (2*a);
                    fb = powl(b,2) * cosl(b) - (2*b);
                    aux = p;
                    p = (a+b)/2;
                    fp = powl(p,2) * cosl(p) - (2*p);
                    Er = fabsl((p - aux)/p);
                    break;
                case 1:
                    fa = (6 - (2/powl(a,2))) * (expl(2 + a)/4) + 1;
                    fb = (6 - (2/powl(b,2))) * (expl(2 + b)/4) + 1;
                    aux = p;
                    p = (a+b)/2;
                    fp = (6 - (2/powl(p,2))) * (expl(2 + p)/4) + 1;
                    Er = fabsl((p - aux)/p);
                    break;
                case 2:
                    fa = powl(a,3) - 3*sinl(powl(a,2)) + 1;
                    fb = powl(b,3) - 3*sinl(powl(b,2)) + 1;
                    aux = p;
                    p = (a+b)/2;
                    fp = powl(p,3) - 3*sinl(powl(p,2)) + 1;
                    Er = fabsl((p - aux)/p);
                    break;
                case 3:
                    fa = powl(a,3) + 6*powl(a,2) + 9.4*a + 2.5;
                    fb = powl(b,3) + 6*powl(b,2) + 9.4*b + 2.5;
                    aux = p;
                    p = (a+b)/2;
                    fp = powl(p,3) + 6*powl(p,2) + 9.4*p + 2.5;
                    Er = fabsl((p - aux)/p);
                    break;
            }

            //Llena la tabla con la informacion de cada iteracion
            vector<string> row;
            row.push_back(to_string(i));
            row.push_back(to_string(a));
            row.push_back(to_string(b));
            row.push_back(to_string(fa));
            row.push_back(to_string(fb));
            row.push_back(to_string(p));
            row.push_back(to_string(fp));
            row.push_back(to_string(Er));
            data.push_back(row);

            //Verifica el criterio
            if(Er < tol){
                this->raiz = to_string(p);
                this->iteracion = to_string(i);
                break;
            }else{
                this->raiz = "No se ha encontrado aún";
                this->iteracion = "-";
            }

            if(fa * fp < 0){
                a = a;
            }else {
                a = p;
            }
            if(fb * fp < 0){
                b = b;
            }else {
                b = p;
            }
        }
    }

    void llenarTablaSecante(unsigned int funcion, double k0, double k, double tol, int iter){
        //Limpiar la tabla
        data.clear();
        vector<string> headers = {"i", "k-1", "k", "f(k-1)", "f(k)", "k+1", "f(k+1)", "Er"};
        data.push_back(headers);

        //Funciones auxiliares
        long double fk0 = 0, fk = 0, k1 = 0, fk1 = 0, Er = 0, aux = 0;

        //Llenar la tabla
        for (int i = 1; i < iter ; ++i) {
            //Verifica la funcion
            switch(funcion){
                case 0:
                    fk0 = powl(k0,2) * cosl(k0) - (2*k0);
                    fk = powl(k,2) * cosl(k) - (2*k);
                    aux = k1;
                    k1 = k - fk * ((k - k0)/(fk - fk0));
                    fk1 = powl(k1,2) * cosl(k1) - (2*k1);
                    Er = fabsl((k1 - aux)/k1);
                    break;
                case 1:
                    fk0 = (6 - (2/powl(k0,2))) * (expl(2 + k0)/4) + 1;
                    fk = (6 - (2/powl(k,2))) * (expl(2 + k)/4) + 1;
                    aux = k1;
                    k1 = k - fk * ((k - k0)/(fk - fk0));
                    fk1 = (6 - (2/powl(k1,2))) * (expl(2 + k1)/4) + 1;
                    Er = fabsl((k1 - aux)/k1);
                    break;
                case 2:
                    fk0 = powl(k0,3) - 3*sinl(powl(k0,2)) + 1;
                    fk = powl(k,3) - 3*sinl(powl(k,2)) + 1;
                    aux = k1;
                    k1 = k - fk * ((k - k0)/(fk - fk0));
                    fk1 = powl(k1,3) - 3*sinl(powl(k1,2)) + 1;
                    Er = fabsl((k1 - aux)/k1);
                    break;
                case 3:
                    fk0 = powl(k0,3) + 6*powl(k0,2) + 9.4*k0 + 2.5;
                    fk = powl(k,3) + 6*powl(k,2) + 9.4*k + 2.5;
                    aux = k1;
                    k1 = k - fk * ((k - k0)/(fk - fk0));
                    fk1 = powl(k1,3) + 6*powl(k1,2) + 9.4*k1 + 2.5;
                    Er = fabsl((k1 - aux)/k1);
                    break;
            }

            //Llena la tabla con la informacion de cada iteracion
            vector<string> row;
            row.push_back(to_string(i));
            row.push_back(to_string(k0));
            row.push_back(to_string(k));
            row.push_back(to_string(fk0));
            row.push_back(to_string(fk));
            row.push_back(to_string(k1));
            row.push_back(to_string(fk1));
            row.push_back(to_string(Er));
            data.push_back(row);

            //Verifica el criterio
            if(Er < tol){
                this->raiz = to_string(k1);
                this->iteracion = to_string(i);
                break;
            }else{
                this->raiz = "No se ha encontrado aún";
                this->iteracion = "-";
            }

            k = k1;
        }
    }

    // Método para mostrar la tabla
    Element Render() {
        Elements table_elements;
        
        // Construir cada fila
        for (const auto& row : data) {
            Elements row_elements;
            for (const auto& cell : row) {
                // Hacer que cada celda se expanda horizontalmente para ocupar todo el ancho
                row_elements.push_back(text(cell) | size(WIDTH, EQUAL, 10) | border | xflex);
            }
            table_elements.push_back(hbox(row_elements));
        }

        // Crear la tabla con bordes y devolver el elemento
        return vbox(table_elements)| size(WIDTH, EQUAL, 150);
    }
};

int main(){

    auto screen = ScreenInteractive::TerminalOutput();

    // Opciones del programa
    /* Radio box */
    int metodo = 0;
    vector<string> metodos = {"Bisección","Secante"};
    auto seleccionMetodos = Radiobox(&metodos, &metodo);

    int funcion = 0;
    vector<string> funciones = {
        "f(x) = x^2 cos(x) - 2x ",
        "f(x) = (6 - (2/x^2))(e^{2 + x}/ 4) + 1",
        "f(x) = x^3 - 3sen(x^2) + 1",
        "f(x) = x^3 + 6x^2 + 9.4x + 2.5",
    };
    auto seleccionFunciones = Radiobox(&funciones, &funcion);
    
    /* Inputs */
    string intervaloInicial;
    string intervaloFinal;
    string tolerancia;
    string numeroIteraciones;
    Component input_intervalo_inicial = Input(&intervaloInicial, "a");
    Component input_intervalo_final = Input(&intervaloFinal, "b");
    Component input_tolerancia = Input(&tolerancia, "E");
    Component input_iteraciones = Input(&numeroIteraciones, "i");

    /* Botón */
    Tabla resultado(5,8);
    auto generarTabla = Button("Generar tabla", [&]{
        switch(metodo){
            case 0:
                resultado.llenarTablaBiseccion(funcion, stod(intervaloInicial), stod(intervaloFinal), stod(tolerancia), stoi(numeroIteraciones));
                break;
            case 1:
                resultado.llenarTablaSecante(funcion, stod(intervaloInicial), stod(intervaloFinal), stod(tolerancia), stoi(numeroIteraciones));
                break;
        }
    }, ButtonOption::Animated(Color::CornflowerBlue));
    
    //Elementos interactivos
    auto components = Container::Vertical({
        seleccionMetodos,
        seleccionFunciones,
        input_intervalo_inicial,
        input_intervalo_final,
        input_tolerancia,
        input_iteraciones,
        generarTabla
    });

    //Render
    auto renderer = Renderer(components, [&] { 
        return window(text("Programa Uno") | center | bold, {
            vbox({
                separatorEmpty(),
                hbox({
                    window(text("Seleccione el método") | center,{
                        seleccionMetodos->Render() | center,
                    }) | size(WIDTH, EQUAL, 40),
                    window(text("Seleccione la función") | center,{
                        seleccionFunciones->Render() | center,
                    }) | size(WIDTH, EQUAL, 60),
                }) | center,
                flexbox({
                    window(text("Intervalo inicial") | center, {
                        input_intervalo_inicial->Render() | size(WIDTH, LESS_THAN, 20) | size(HEIGHT, LESS_THAN, 3)
                    }),
                    window(text("Intervalo final") | center, {
                        input_intervalo_final->Render() | size(WIDTH, LESS_THAN, 20) | size(HEIGHT, LESS_THAN, 3)
                    }),
                    window(text("Tolerancia") | center, {
                        input_tolerancia->Render() | size(WIDTH, LESS_THAN, 20) | size(HEIGHT, LESS_THAN, 3)
                    }),
                    window(text("Número de iteraciones") | center, {
                        input_iteraciones->Render() | size(WIDTH, LESS_THAN, 20) | size(HEIGHT, LESS_THAN, 3)
                    })
                }, FlexboxConfig()
                    .Set(FlexboxConfig::JustifyContent::SpaceEvenly)
                ),
                generarTabla->Render() | center,
                separatorEmpty(),
                window(text("Tabla de resultados") | center | color(Color::Gold1),{
                    resultado.Render() | center,
                }) | center,
                separatorEmpty(),
                window(text("Resultados") | center | color(Color::Gold1),{
                    text("Raíz: " + resultado.raiz + " en el intervalo [" + intervaloInicial + "," + intervaloFinal + "] en la iteración: " + resultado.iteracion) | center | color(Color::Gold1)
                })| center,
            })
        }) | size(WIDTH, LESS_THAN, 150) | center;
    });

    screen.Loop(renderer);

    return 0;
}