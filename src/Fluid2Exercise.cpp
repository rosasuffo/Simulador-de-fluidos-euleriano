#include "Scene.h"

#include "Numeric/PCGSolver.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace asa
{
namespace
{
Vector3 bilinearInt_Vector3(Array2<Vector3> values, Vector2 idx)
{
    int x0(floor(idx.x));
    int x1(ceil(idx.x));
    int y0(floor(idx.y));
    int y1(ceil(idx.y));
    float deltaX = (x0 != x1) ? idx.x - x0 : 0.f;
    float deltaY = (y0 != y1) ? idx.y - y0 : 0.f;

    Vector3 f_x0y0 = values[Index2(x0, y0)];
    Vector3 f_x0y1 = values[Index2(x0, y1)];
    Vector3 f_x1y0 = values[Index2(x1, y0)];
    Vector3 f_x1y1 = values[Index2(x1, y1)];

    Vector3 q0 = f_x0y0 * (1.f - deltaX) + f_x1y0 * deltaX;
    Vector3 q1 = f_x0y1 * (1.f - deltaX) + f_x1y1 * deltaX;
    return q0 * (1.f - deltaY) + q1 * deltaY;
}

float bilinearInt_float(Array2<float> values, Vector2 idx)
{
    int x0(floor(idx.x));
    int x1(ceil(idx.x));
    int y0(floor(idx.y));
    int y1(ceil(idx.y));
    float deltaX = (x0 != x1) ? idx.x - x0 : 0.f;
    float deltaY = (y0 != y1) ? idx.y - y0 : 0.f;

    float f_x0y0 = values[Index2(x0, y0)];
    float f_x0y1 = values[Index2(x0, y1)];
    float f_x1y0 = values[Index2(x1, y0)];
    float f_x1y1 = values[Index2(x1, y1)];

    float q0 = f_x0y0 * (1.f - deltaX) + f_x1y0 * deltaX;
    float q1 = f_x0y1 * (1.f - deltaX) + f_x1y1 * deltaX;
    return q0 * (1.f - deltaY) + q1 * deltaY;
}

////////////////////////////////////////////////
// Add any reusable classes or functions HERE //
////////////////////////////////////////////////
}  // namespace


// advection
void Fluid2::fluidAdvection(const float dt)
{
    // La idea del algoritmo de adveccion
    // 1. Calcular posicion donde se situa la propiedad a calcular (velocidad o tinta) en el paso anterior. 
    // Para ello, calcular la velocidad en el punto actual y aplicarla de manera inversa a la posicion.
    // 2. Interpolar valor de la propiedad en la posicion obtenida del paso anterior.
    // 3. Valor interpolado de la propiedad = resultado que buscamos.

    Array2<float> velXCopy = this->velocityX;
    Array2<float> velYCopy = this->velocityY;

    // Ink advecion HERE
    Array2<Vector3> inkRGBCopy = this->inkRGB;
    for (uint i = 0; i < Scene::nCellsX; i++) {
        for (uint j = 0; j < Scene::nCellsY; j++) {
            Index2 ij(i, j);
            Vector2 x = this->grid.getCellPos(ij);

            // se obtiene el vector velocidad en el punto medio. para ello se interpola el de cada cara de la rejilla
            float v_x = (velXCopy[Index2(i, j)] + velXCopy[Index2(i + 1, j)]) * 0.5f;
            float v_y = (velYCopy[Index2(i, j)] + velYCopy[Index2(i, j + 1)]) * 0.5f;
            Vector2 v(v_x, v_y);

            // ecuacion de cinematica para interpolacion de posiciones, con la velocidad invertida para saber el paso anterior
            Vector2 xPrev = x - v * dt;
            Vector2 ijPrev = this->grid.getCellIndex(xPrev);

            // se clampean las posiciones para que no se salgan del dominio de la rejilla
            ijPrev.x = std::max(0.f, std::min(ijPrev.x, (float)Scene::nCellsX - 1.001f));
            ijPrev.y = std::max(0.f, std::min(ijPrev.y, (float)Scene::nCellsY - 1.001f));

            Vector3 ink = bilinearInt_Vector3(inkRGBCopy, ijPrev);
            this->inkRGB[Index2(i, j)] = ink;
            }
    }

    // Velocity advection term HERE
    Index2 facesX_size = this->grid.getSizeFacesX();
    for (uint i = 0; i < facesX_size.x; i++) {
        for (uint j = 0; j < facesX_size.y; j++) {
            Index2 ij(i, j);
            Vector2 x = this->grid.getFacePosX(ij);

            // se busca el valor de la velocidad en el punto medio de cada cara izquierda de la rejilla
            // ese valor de la velocidad x ya lo tenemos porque es ahi donde se guarda. pero la velocidad en y
            // se guarda en las caras inferiores de la rejilla. hay que interpolar con los valores de la
            // velocidad en Y de: la celda actual, la celda izda, la celda superior y la celda izda superior.
            float vx        = velXCopy[Index2(i, j)];
            float vy_left   = (i == 0) ? 0.f : velYCopy[Index2(i - 1, j)];
            float vy_center = (i == facesX_size.x - 1) ? 0.f : velYCopy[Index2(i, j)]; // vy tiene facesX_size.x - 1 elementos
            float vy_leftup = (i == 0) ? 0.f : velYCopy[Index2(i - 1, j + 1)];
            float vy_up     = (i == facesX_size.x -1) ? 0.f : velYCopy[Index2(i, j + 1)];
            float vy        = (vy_left + vy_center + vy_leftup + vy_up) * 0.25f;
            Vector2 v(vx, vy);

            Vector2 xPrev = x - v * dt;
            Vector2 ijPrev = this->grid.getFaceIndexX(xPrev);

            ijPrev.x = std::max(0.f, std::min(ijPrev.x, (float)facesX_size.x - 1.001f));
            ijPrev.y = std::max(0.f, std::min(ijPrev.y, (float)facesX_size.y - 1.001f));

            float new_vx = bilinearInt_float(velXCopy, ijPrev);
            this->velocityX[Index2(i, j)] = new_vx;
        }
    }

    Index2 facesY_size = this->grid.getSizeFacesY();
    for (uint i = 0; i < facesY_size.x; i++) {
        for (uint j = 0; j < facesY_size.y; j++) {
            Index2 ij(i, j);
            Vector2 x = this->grid.getFacePosY(ij);

            float vx_down       = (j == 0) ? 0.f : velXCopy[Index2(i, j - 1)];
            float vx_center     = (j == facesY_size.y - 1) ? 0.f : velXCopy[Index2(i, j)];
            float vx_downright  = (j == 0) ? 0.f : velXCopy[Index2(i + 1, j - 1)];
            float vx_right      = (j == facesY_size.y - 1) ? 0.f : velXCopy[Index2(i + 1, j)];
            float vx            = (vx_down + vx_center + vx_downright + vx_right) * 0.25f;
            float vy            = velYCopy[Index2(i, j)];
            Vector2 v(vx, vy);

            Vector2 xPrev = x - v * dt;
            Vector2 ijPrev = this->grid.getFaceIndexY(xPrev);

            ijPrev.x = std::max(0.f, std::min(ijPrev.x, (float)facesY_size.x - 1.001f));
            ijPrev.y = std::max(0.f, std::min(ijPrev.y, (float)facesY_size.y - 1.001f));

            float velY = bilinearInt_float(velYCopy, ijPrev);
            this->velocityY[Index2(i, j)] = velY;
        }
    }
}

void Fluid2::fluidEmission()
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Emitters contribution HERE
        // Definimos unas formas en nuestro dominio (por ejemplo, cajas).
        // Toda celdilla que caiga en esas cajas que son emisores hay que ponerles ink y vel concretas
        // Para el resultado del pdf mirar el documento subido con el dominio, velocidad e ink de cada emisor

        /*float emit_center = 2;
        float emit_width = 1;
        float emit_height = 1;
        const float ink = 1;
        float vel = 1;

        for (int i = emit_center - emit_width * 0.5f; i < emit_center + emit_width * 0.5f; i++) {
            for (int j = emit_center - emit_height * 0.5f; j < emit_center + emit_height; j++) {
                Index2 idx(i, j);
                this->inkRGB[idx] = Vector3(ink, ink, ink);
                this->velocityY[idx] = vel;
            }
        }*/ 
        
        // dimension de cada cuadrado de color = tamaño de la rejilla / 10
        float emit_width = Scene::nCellsX / 10;
        float emit_height = Scene::nCellsY / 30;
        float emit_center = Scene::nCellsX / 2;

        float init_cell_x = emit_center - (emit_width * 1.5);
        
        for (int j = 1; j <= emit_height; j++) {  // para cada fila, de abajo a arriba de la rejilla
            for (int i = init_cell_x; i <= emit_center + (emit_width * 1.5); i++) {
                if (i < init_cell_x + emit_width)
                    this->inkRGB[Index2(i,j)] = Vector3(1, 0, 0);
                else if (i < init_cell_x + emit_width * 2)
                    this->inkRGB[Index2(i, j)] = Vector3(0, 1, 0);
                else
                    this->inkRGB[Index2(i, j)] = Vector3(0, 0, 1);

                this->velocityY[Index2(i,j)] = 8;
            }
        }  
    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Gravity term HERE
        // Como es la gravedad, �nicamente se aplica a velocidades en y
        // Foreach velocityY i, j
        // Resolver expl�citamente con g = Scene::kGravity
        // Si nos ponemos creativos, se puede hacer similar a los emisores: a�adir dominios de viento...

        Index2 faceSizeY = grid.getSizeFacesY();

        for (uint i = 0; i < faceSizeY.x; i++) {
            for (uint j = 0; j < faceSizeY.y; j++) {
                this->velocityY[Index2(i, j)] += Scene::kGravity * dt;
            }
        }
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    // Viscosity term HERE
    
    // se resuelve la formula del termino viscoso de navier-strokes con una aproximacion explicita por diferencias finitas
    // v' = v + (dt / densidad) * viscosidad * ((v_right - 2*v + v_left) / dx^2 + (v_up - 2*v + v_down) / dy^2)
    // hay que descomponerla para la velocidad en x y la velocidad en y.
    // Condiciones de contorno -> se clampea a la parte conocida (por ejemplo, v[-1] -> v[0])
    
    // una copia de la velocidad para actualizar y otra de lectura
    Array2<float> velX_write = this->velocityX;
    Array2<float> velY_write = this->velocityY;        
    Array2<float> velX_read = velX_write;
    Array2<float> velY_read = velY_write;
    uint Nx = Scene::nCellsX;
    uint Ny = Scene::nCellsY;
    float deltaX = this->getGrid().getDx().x;
    float deltaY = this->getGrid().getDx().y;
    float kX = (dt * Scene::kViscosity / Scene::kDensity) / (deltaX * deltaX);
    float kY = (dt * Scene::kViscosity / Scene::kDensity) / (deltaY * deltaY);
    
    Index2 facesX_size = this->grid.getSizeFacesX();
    Index2 facesY_size = this->grid.getSizeFacesY();
    
    for (uint i = 1; i < Nx; i++) {
        for (uint j = 0; j < Ny; j++) {
            float vx = velX_read[Index2(i, j)];
            // left
            float vx_left = velX_read[Index2(i - 1, j)];
            // right
            float vx_right = velX_read[Index2(i + 1, j)];
            // down
            float vx_down = (j > 0) ? velX_read[Index2(i, j - 1)] : vx;
            // up
            float vx_up = (j < Ny - 1) ? velX_read[Index2(i, j + 1)] : vx;
    
            float vx_new = kX * (vx_right - 2 * vx + vx_left) + kY * (vx_up - 2 * vx + vx_down);
            
            velX_write[Index2(i, j)] += vx_new;
        }
    }
    
    for (uint i = 0; i < Nx; i++) {
        for (uint j = 1; j < Ny; j++) {
            float vy = velY_read[Index2(i, j)];
            // left
            float vy_left = (i > 0) ? velY_read[Index2(i - 1, j)] : vy;
            // right
            float vy_right = (i < Nx - 1) ? velY_read[Index2(i + 1, j)] : vy;
            // down
            float vy_down = velY_read[Index2(i, j - 1)];
            // right
            float vy_up = velY_read[Index2(i, j + 1)];
    
            float vy_new = kX * (vy_right - 2 * vy + vy_left) + kY * (vy_up - 2 * vy + vy_down);
    
            velY_write[Index2(i, j)] += vy_new;
        }
    }
    
    this->velocityX = velX_write;
    this->velocityY = velY_write;
}

void Fluid2::fluidPressureProjection(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {

        // Incompressibility / Pressure term HERE
        int Nx = Scene::nCellsX;
        int Ny = Scene::nCellsY;
        Array2<float> vx = this->velocityX;
        Array2<float> vy = this->velocityY;
        // se asume que dx = dy
        float dx = this->getGrid().getDx().x;
        float K = -Scene::kDensity / (dt * dx);

        // Condiciones de contorno: Set velocity components in all boundaries to 0
        // Foreach Y (row) make vX0 (first column), vXN(last column) = 0
        // Foreach X make vY0, vYN = 0
        for (uint i = 0; i < Nx; i++) { // para cada elemento de una fila (nx+1)
            vy[Index2(i, 0)] = 0;  // primera fila
            vy[Index2(i, Ny)] = 0; // ult fila
        }
        for (uint i = 0; i < Ny; i++) { // para cada elemento de una columna (Ny+1)
            vx[Index2(0, i)] = 0;   // primera columna
            vx[Index2(Nx, i)] = 0;  // ult columna
        }

        // Inicializar el solver
        // Elegir precision y tolerancia dependiendo de si queremos velocidad o exactitud
        // float => tolerance_factor = 1e-3, iterations = 200   
        // double => tolerance_factor = 1e-6, iterations = 200
        PCGSolver<double> solver;
        solver.set_solver_parameters(1e-2, 1000);

        // SISTEMA DE ECUACIONES A RESOLVER: 
        
        // RHS (vector b del sistema de ecuaciones) -> divergencia del campo de velocidades
        std::vector<double> rhs(Nx * Ny, 0.0);
        // Matriz A -> laplaciano de la presion
        SparseMatrix<double> A(Nx * Ny, 5);
        // incognita = presiones
        std::vector<double> p(Nx * Ny, 0.0);

        for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
        p[i * Ny + j] = (double)this->pressure[Index2(i, j)];
    }
}

        // Rellenar matriz A y RHS
        float dxdx = dx * dx;
        A.zero();
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                int idx = i * Ny + j;

                // RHS =  - (dens / dt) * ((v_up - v) / dx + (v_left - v) / dy)
                // Asumiendo que dx = dy y separando la ecuación en dos partes, K * Divergence(v_ij)
                // Divergence(v_ij) = (v_up - v) + (v_left - v)
                // K = dens / (dt * dx);
                float div = (vx[Index2(i + 1, j)] - vx[Index2(i, j)]) + (vy[Index2(i, j + 1)] - vy[Index2(i, j)]);
                rhs[idx] = (double)(K * div);

                // matriz A: operador laplaciano negativo 
                // Asumiendo que dx = dy, queda: A = (p_left + p_right + p_down + p_up - 4*p) / dx^2
                // Si aislamos los valores de p de la ecuación, nos queda una matriz que hay que rellenar de la siguiente manera.
                // vecino arriba = -1 (si existe)
                // vecino abajo = -1 (si existe)
                // vecino dcha = -1 (si existe)
                // vecino izda = -1 (si existe)
                // elemento central = nº de vecinos 
                // use add_to_element, set_element

                int neighbors = 0;
                if (i > 0) { // vecino abajo
                    A.set_element(idx, (i - 1) * Ny + j, -1.0/ (dx * dx));
                    neighbors++;
                }
                if (i < Nx - 1) { // vecino arriba
                    A.set_element(idx, (i + 1) * Ny + j, -1.0 / (dx * dx));
                    neighbors++;
                }
                if (j > 0) { // vecino izda
                    A.set_element(idx, i * Ny + (j - 1), -1.0 / (dx * dx));
                    neighbors++;
                }
                if (j < Ny - 1) { // vecino dcha
                    A.set_element(idx, i * Ny + (j + 1), -1.0 / (dx * dx));
                    neighbors++;
                }

                A.set_element(idx, idx, (double)neighbors / (dx * dx));
            }
        }

        // Solve
        // Se pueden mirar nIter y más datos con T residual_out y iterations_out como parámetros adicionales
        // residual_out: mide la precisión del resultado. si la función devuelve true, significa que residual_out tiene un valor menor que una tolerancia establecida
        // iterations_out: nº de iteraciones realizadas para alcanzar el resultado.
        double residual_out = 0;
        int iters_out = 0;
        solver.solve(A, rhs, p, residual_out, iters_out);

        std::cout << "Iteraciones realizadas: " << iters_out << std::endl;

        // Se debe calcular la media de las presiones, para asegurar que las presiones estén siempre centradas en 0 y estabilizar el sistema
        double sum_p = 0.0;
        for (int i = 0; i < Nx * Ny; i++) {
            sum_p += p[i];
        }
        double mean_p = sum_p / (double)(Nx * Ny);

        // Aplicar las P a la rejilla -> pasar todos los valores a Array2 pressure (definido en Fluid2.h)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                this->pressure[Index2(i, j)] = (float)(p[i * Ny + j] - mean_p);

        // Aplicar últimas fórmulas del pdf -> actualizar todas las velocidades del sistema con las presiones
        float p_K = dt / (Scene::kDensity * dx);
        // las velocidades en los extremos valen 0, además en vx, sus últimas columnas no tienen acceso a valores de presion porque no existen
        //                                                    vy, sus últimas filas no tienen acceso a valores de presion porque no existen
        
        //  -----------
        // |-> |-> |-> |->  2           para vx, se recorre la rejilla desde la columna 
        //  -----------                 1 a la 2 y desde la fila 0 a la 2
        // |-> |-> |-> |->  1           
        //  -----------                 
        // |-> |-> |-> |->  0
        //  -----------
        //   0   1   2   3   
        for (int i = 1; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                vx[Index2(i, j)] -= p_K * (this->pressure[Index2(i, j)] - this->pressure[Index2(i-1, j)]);

        //   ^   ^   ^
        //   |   |   |    3
        //  -----------                 para vy, se recorre la rejilla desde la columna
        // | ^ | ^ | ^ |                0 a la 2 y desde la fila 1 a la 2
        // | | | | | | |  2      
        //  -----------            
        // | ^ | ^ | ^ |
        // | | | | | | |  1
        //  -----------
        // | ^ | ^ | ^ |
        // | | | | | | |  0
        //  -----------
        //   0   1   2   
        for (int i = 0; i < Nx; i++)
            for (int j = 1; j < Ny; j++)
                vy[Index2(i, j)] -= p_K * (this->pressure[Index2(i, j)] - this->pressure[Index2(i, j-1)]);

        // Se puede comprobar si para cada una de las celdas la fórmula nabla·uij = 0, en cuyo caso estaría bien hecho

        this->velocityX = vx;
        this->velocityY = vy;
    }
}
}  // namespace asa
