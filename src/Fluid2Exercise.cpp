#include "Scene.h"

#include "Numeric/PCGSolver.h"

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
    float deltaY = (y0 != y1) ? idx.y - y0: 0.f;

    Vector3 f_x0y0 = values[x0, y0];
    Vector3 f_x0y1 = values[x0, y1];
    Vector3 f_x1y0 = values[x1, y0];
    Vector3 f_x1y1 = values[x1, y1];

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

    float f_x0y0 = values[x0, y0];
    float f_x0y1 = values[x0, y1];
    float f_x1y0 = values[x1, y0];
    float f_x1y1 = values[x1, y1];

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
    {
        // Ink advecion HERE
        Array2<Vector3> inkRGBCopy = this->inkRGB;
        Array2<float> velXCopy = this->velocityX;
        Array2<float> velYCopy = this->velocityY;

        for (uint i = 0; i < Scene::nCellsX; i++) {
            for (uint j = 0; j < Scene::nCellsY; j++) {
                Index2 ij(i, j);
                Vector2 x = this->grid.getCellPos(ij);
                //Vector2 v = Vector2(velXCopy[i, j], velYCopy[i, j]);
                float u1 = (velXCopy[i, j] + velXCopy[i + 1, j]) * 0.5f;
                float u2 = (velYCopy[i, j] + velYCopy[i, j + 1]) * 0.5f;
                Vector2 v(u1, u2);

                Vector2 xPrev = x - v * dt;
                Vector2 ijPrev = this->grid.getCellIndex(xPrev);

                ijPrev.x = std::max(0.f, std::min(ijPrev.x, (float)Scene::nCellsX - 1));
                ijPrev.y = std::max(0.f, std::min(ijPrev.y, (float)Scene::nCellsY- 1)) ;

                Vector3 ink = bilinearInt_Vector3(inkRGBCopy,ijPrev);
                this->inkRGB[i, j] = ink;
            }
        }
    }
    {
        // Velocity advection term HERE
        Array2<float> velXCopy = this->velocityX;
        Array2<float> velYCopy = this->velocityY;
        Array2<float> velXNew = velXCopy;
        Array2<float> velYNew = velYCopy;

        Index2 facesX_size = this->grid.getSizeFacesX();

        for (uint i = 0; i < facesX_size.x; i++) {
            for (uint j = 0; j < facesX_size.y; j++) {
                Index2 ij(i, j);
                Vector2 x = this->grid.getFacePosX(ij);

                float u1 = velXCopy[i, j];
                float left = (i == 0) ? 0.f : velYCopy[i - 1, j];
                float leftup = (i == 0) ? 0.f : velYCopy[i - 1, j + 1];
                float u2 = (left + velYCopy[i, j] + leftup + velYCopy[i, j + 1]) * 0.25f;
                Vector2 v(u1, u2);

                Vector2 xPrev = x - v * dt;
                Vector2 ijPrev = this->grid.getFaceIndexX(xPrev);

                ijPrev.x = std::max(0.f, std::min(ijPrev.x, (float)Scene::nCellsX));
                ijPrev.y = std::max(0.f, std::min(ijPrev.y, (float)Scene::nCellsY - 1));

                float velX = bilinearInt_float(velXCopy, ijPrev);
                velXNew[i, j] = velX;  // Escritura segura
            }
        }

        Index2 facesY_size = this->grid.getSizeFacesY();

        for (uint i = 0; i < facesY_size.x; i++) {
            for (uint j = 0; j < facesY_size.y; j++) {  // Límite corregido
                Index2 ij(i, j);
                Vector2 x = this->grid.getFacePosY(ij);

                float u1 = velYCopy[i, j];
                float down = (j == 0) ? 0.f : velXCopy[i, j - 1];
                float downright = (j == 0) ? 0.f : velXCopy[i + 1, j - 1];
                float u2 = (down + velYCopy[i, j] + downright + velYCopy[i + 1, j]) * 0.25f;
                Vector2 v(u1, u2);

                Vector2 xPrev = x - v * dt;
                Vector2 ijPrev = this->grid.getFaceIndexY(xPrev);

                ijPrev.x = std::max(0.f, std::min(ijPrev.x, (float)Scene::nCellsX - 1));
                ijPrev.y = std::max(0.f, std::min(ijPrev.y, (float)Scene::nCellsY));

                float velY = bilinearInt_float(velYCopy, ijPrev);
                velYNew[i, j] = velY;  // Escritura segura
            }
        }

        this->velocityX = velXNew;
        this->velocityY = velYNew;
    }
}

void Fluid2::fluidEmission()
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Emitters contribution HERE
        // Definimos unas formas en nuestro dominio (por ejemplo, cajas).
        // Toda celdilla que caiga en esas cajas que son emisores hay que ponerles ink y vel concretas
        // Para el resultado del pdf mirar el documento subido con el dominio, velocidad e ink de cada emisor

        float emit_center = 2;
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
        }
    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Gravity term HERE
        // Como es la gravedad, únicamente se aplica a velocidades en y
        // Foreach velocityY i, j
        // Resolver explícitamente con g = Scene::kGravity
        // Si nos ponemos creativos, se puede hacer similar a los emisores: añadir dominios de viento...

        for (uint i = 0; i < Scene::nCellsX; i++) {
            for (uint j = 0; j < Scene::nCellsY; j++) {
                this->velocityY[i, j] += Scene::kGravity * dt;
            }
        }
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Viscosity term HERE
        // La ecuación del pdf son en realidad 2: una para x y otra para y
        // Copiar velX y velY
        // Aplicar término viscoso a velX y velY:
        // Foreach face i, j
        // Diferencias finitas centradas de segundo orden (fórmula pdf); hay métodos en grid??
        // Aplicar integración de tiempo
        // Multiplicar por término viscoso y dividir por densidad
        // Condiciones de contorno -> se clampea a la parte conocida (por ejemplo, v[-1] -> v[0])

        Array2<float> velX = this->velocityX;
        Array2<float> velY = this->velocityY;        
        Array2<float> velXOld = velX;
        Array2<float> velYOld = velY;
        uint Nx = Scene::nCellsX;
        uint Ny = Scene::nCellsY;
        float deltaX = this->getGrid().getDx().x;
        float deltaY = this->getGrid().getDx().y;
        float kX = (dt * Scene::kViscosity / Scene::kDensity) / (deltaX * deltaX);
        float kY = (dt * Scene::kViscosity / Scene::kDensity) / (deltaY * deltaY);

        for (uint i = 0; i <= Nx; i++) {
            for (uint j = 0; j <= Ny; j++) {
                float vel = velXOld[i, j];
                // left
                Index2 idx_left = (i > 0) ? Index2(i - 1,j) : Index2(i, j);
                float vel_left = velXOld[idx_left];
                // right
                Index2 idx_right = (i <= Nx) ? Index2(i + 1,j) : Index2(i, j);
                float vel_right = velXOld[idx_right];
                // up
                Index2 idx_down = (j > 0) ?  Index2(i, j-1) :  Index2(i, j);
                float vel_down = velXOld[idx_down];
                // down
                Index2 idx_up = (j <= Ny) ? Index2(i, j + 1) : Index2(i, j);
                float vel_up = velXOld[idx_up];

                float new_vel = kX * (vel_right - 2 * vel + vel_left) + kY * (vel_up - 2 * vel + vel_down);
                
                if (i >= 0 && i < Nx - 1 && j > 0 && j < Ny) {
                    velX[i, j] += new_vel;
                }
            }
        }

        for (uint i = 0; i <= Nx; i++) {
            for (uint j = 0; j <= Ny; j++) {
                float vel = velYOld[i, j];
                // left
                Index2 idx_left = (i > 0) ? Index2(i - 1, j) : Index2(i, j);
                float vel_left = velYOld[idx_left];
                // right
                Index2 idx_right = (i < Nx - 1) ? Index2(i + 1, j) : Index2(i, j);
                float vel_right = velYOld[idx_right];
                // up
                Index2 idx_down = (j > 0) ? Index2(i, j - 1) : Index2(i, j);
                float vel_down = velYOld[idx_down];
                // down
                Index2 idx_up = (j < Ny) ? Index2(i, j + 1) : Index2(i, j);
                float vel_up = velYOld[idx_up];

                float new_vel = kX * (vel_right - 2 * vel + vel_left) + kY * (vel_up - 2 * vel + vel_down);

                if (i >= 0 && i < Nx && j > 0 && j < Ny) {
                    velY[i, j] += new_vel;
                }
            }
        }

        this->velocityX = velX;
        this->velocityY = velY;
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        int Nx = Scene::nCellsX;
        int Ny = Scene::nCellsY;
        Array2<float> velocityX = this->velocityX;
        Array2<float> velocityY = this->velocityY;
        // se asume que dx = dy
        float dx = this->getGrid().getDx().x;
        //float K = -1.f / (dt * dx);
        float K = -(dx * Scene::kDensity) / dt;

        // Incompressibility / Pressure term HERE
        // Set normal velocity components in all boundaries to 0
        // Foreach Y make vX0, vXN = 0
        // Foreach X make vY0, vYN = 0
        for (uint i = 0; i < Nx; i++) {
            velocityY[i, 0] = 0;
            velocityY[i, Ny] = 0;
        }
        for (uint i = 0; i < Ny; i++) {
            velocityX[0, i] = 0;
            velocityX[Nx, i] = 0;
        }

        // En solver elegir precisión y tolerancia dependiendo de si queremos velocidad o exactitud
        // float => tolerance_factor = 1e-3, iterations = 200
        // double => tolerance_factor = 1e-6, iterations = 200
        // PCGSolver solver;
        // solver.set_solver_parameters()
        PCGSolver<double> solver;
        solver.set_solver_parameters(1e-6, 200);

        // RHS (vector b del sistema de ecuaciones) -> divergencia del campo de velocidades
        // Foreach i, j Apply K * Divergence(Uij)
        std::vector<double> rhs(Nx * Ny, 0.f);

        // Matriz A -> laplaciano de la presion
        // Foreach i, j
        //   Apply laplacian(P) -> Ha dicho algo sobre traducir de i j a no sé qué, tener cuidado con índices
        //   use add_to_element, set_element
        // Mirar asunción delta x = delta y; en el enumerador se queda -4pij, aunque se resta 1 por contorno adyacente
        SparseMatrix<double> A(Nx * Ny);

        // P
        std::vector<double> p(Nx * Ny, 0.f);

        // Rellenar matriz A y RHS
        A.zero(); 

        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = i + j * Nx;

                float div = (velocityX[i + 1, j] - velocityX[i, j]) + (velocityY[i, j + 1] - velocityY[i, j]);
                rhs[idx] = (double)(K * div);

                int neighbors = 0;
                if (i > 0) {
                    A.set_element(idx, (i - 1) + j * Nx, -1.0);
                    neighbors++;
                }
                if (i < Nx - 1) { 
                    A.set_element(idx, (i + 1) + j * Nx, -1.0);
                    neighbors++;
                }
                if (j > 0) {
                    A.set_element(idx, i + (j - 1) * Nx, -1.0);
                    neighbors++;
                }
                if (j < Ny - 1) {
                    A.set_element(idx, i + (j + 1) * Nx, -1.0);
                    neighbors++;
                }

                A.set_element(idx, idx, (double)neighbors);
            }
        }

        // Solve
        // solver.solve(A, RHS, P)
        // Se pueden mirar nIter y más datos con T residual_out y iterations_out como parámetros adicionales
        double residual_out = 0;
        int iters_out = 0;
        solver.solve(A, rhs, p, residual_out, iters_out);

        // Aplicar las P a la rejilla -> pasar todos los valores a Array2 pressure (definido en Fluid2.h)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                this->pressure[i, j] = (float)p[i + j * Nx];

        // Aplicar últimas fórmulas del pdf -> actualizar todas las velocidades del sistema con las presiones
        float p_K = dt / (Scene::kDensity * dx);
        for (int i = 1; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                velocityX[i, j] -= p_K * (this->pressure[i, j] - this->pressure[i - 1, j]);
        for (int i = 0; i < Nx; i++)
            for (int j = 1; j < Ny; j++)
                velocityY[i, j] -= p_K * (this->pressure[i, j] - this->pressure[i, j - 1]);
        // Se puede comprobar si para cada una de las celdas la fórmula nabla·uij = 0, en cuyo caso estaría bien hecho
        this->velocityX = velocityX;
        this->velocityY = velocityY;
    }
}
}  // namespace asa
