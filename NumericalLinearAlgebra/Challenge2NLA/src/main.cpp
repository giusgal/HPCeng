#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <cstdlib>

#include "../include/utils.hpp"

using Image = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using SpImage = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// Just a semantic thing
using Vector = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using SpVector = Eigen::SparseMatrix<double, Eigen::RowMajor>;

int main() {
    /*
    TASK1
    Load the image as an Eigen matrix A with size m×n. Each entry in the matrix corresponds
    to a pixel on the screen and takes a value somewhere between 0 (black) and 255 (white).
    Compute the matrix product A^TA and report the euclidean norm of A^TA.
    */

    Image A;
    Utils::loadImage("../data/Einstein.jpg", A);

    Image ATA{A.transpose()*A};
    double ATAnorm = ATA.norm();
    std::cout << "[INFO] norm(A^TA) = " << ATAnorm << std::endl;

    /*
    TASK2
    Solve the eigenvalue problem A^TAx = λx using the proper solver provided by the Eigen
    library. Report the two largest computed singular values of A.
    */
    
    // A^TA is symmetric => SelfAdjointEigenSolver
    Eigen::SelfAdjointEigenSolver<Image> eigensolver(ATA);
    if(eigensolver.info() != Eigen::Success) {
        std::cerr << "[ERROR] Couldn't compute eigenvalues of ATA" << std::endl;
        exit(1);
    }
    Vector ATAeigenvalues = eigensolver.eigenvalues();
    // std::cout << ATAeigenvalues << std::endl;
    
    // The singular values of A are the (ordered) square roots of the eigenvalues of A^TA,
    //  therefore
    double ATAmaxeigen1 = ATAeigenvalues(ATA.rows()-1);
    double ATAmaxeigen2 = ATAeigenvalues(ATA.rows()-2);
    std::cout << "[INFO] maxeigen1(A^TA) = " << ATAmaxeigen1 
        << std::endl;
    std::cout << "[INFO] maxeigen2(A^TA) = " << ATAmaxeigen2 
        << std::endl;
    double AsingularValue1 = std::sqrt(ATAmaxeigen1);
    double AsingularValue2 = std::sqrt(ATAmaxeigen2);
    std::cout << "[INFO] σ1(A) = " << AsingularValue1 << std::endl;
    std::cout << "[INFO] σ2(A) = " << AsingularValue2 << std::endl;


    // TODO: Check tolerances
    /*
    TASK3 (LIS) [execute solveLis.sh after this file]
    Export matrix ATA in the matrix market format and move it to the lis-2.1.6/test
    folder. Using the proper iterative solver available in the LIS library compute the largest
    eigenvalue of ATA up to a tolerance of 10^−8. Report the computed eigenvalue. Is the result
    in agreement with the one obtained in the previous point?

    Result:
     max_eig(A^TA) = 1.608332e+04 (it's the same)
     relative residual = 1.866013e-09
     iterations = 8
    */
    
    // Export the matrix
    if(!Eigen::saveMarket(ATA, "../data/ATA.mtx")) {
        std::cerr << "[ERROR] Couldn't save matrix ATA in the matrix market format" << std::endl;
        exit(1);
    }


    // TODO: See this better (check tolerances)
    /*
    TASK4 (LIS) [execute solveLis.sh after this file]
    Find a shift μ ∈ R yielding an acceleration of the previous eigensolver. Report μ and the
    number of iterations required to achieve a tolerance of 10^−8.
    
    Result:
     μ = 700
     relative residual = 7.757747e-10
     iterations = 7
    */


    // TODO: Is there a difference in computing the norm when considering the full and reduced SVD?
    /*
    TASK5
    Using the SVD module of the Eigen library, perform a singular value decomposition of the
    matrix A. Report the Euclidean norm of the diagonal matrix Σ of the singular values.
    */

    // Eigen::BDCSVD<Image> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::BDCSVD<Image> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Vector AsingularValues = svd.singularValues();
    Image AU = svd.matrixU();
    Image AV = svd.matrixV();
    // std::cout << "[TEST]" << AsingularValues(0) << std::endl;
    // std::cout << "[TEST]" << AsingularValues(1) << std::endl;
    Image Asigma = AsingularValues.asDiagonal();
    double AsigmaNorm = Asigma.norm();
    std::cout << "[INFO] norm(Σ_A) = " << AsigmaNorm << std::endl;


    // TODO: is this correct? (I mean, the way the matrices are computed and transformed
    //  into sparse matrices)
    // TODO: SEE THIS
    /*
    TASK6
    Compute the matrices C and D described in (1) assuming k = 40 and k = 80. Report the
    number of nonzero entries in the matrices C and D
    */
    double k1 = 40;
    SpImage Ck1 = AU.block(0,0,AU.rows(),k1).sparseView();
    SpImage Dk1 = (AV.block(0,0,AV.rows(),k1)*Asigma.block(0,0,k1,k1)).sparseView();

    double k2 = 80;
    SpImage Ck2 = AU.block(0,0,AU.rows(),k2).sparseView();
    SpImage Dk2 = (AV.block(0,0,AV.rows(),k2)*Asigma.block(0,0,k2,k2)).sparseView();

    int Ck1nnz = Ck1.nonZeros();
    int Dk1nnz = Dk1.nonZeros();
    int Ck2nnz = Ck2.nonZeros();
    int Dk2nnz = Dk2.nonZeros();

    std::cout << "[INFO] nnz(Ck1) = " << Ck1nnz << std::endl;
    std::cout << "[INFO] nnz(Dk1) = " << Dk1nnz << std::endl;
    std::cout << "[INFO] nnz(Ck2) = " << Ck2nnz << std::endl;
    std::cout << "[INFO] nnz(Dk2) = " << Dk2nnz << std::endl;
    

    /*
    TASK7
    Compute the compressed images as the matrix product CD^T (again for k = 40 and k = 80).
    Export and upload the resulting images in .png.
    */
    SpImage AcmpK1 = Ck1*Dk1.transpose();
    SpImage AcmpK2 = Ck2*Dk2.transpose();

    Utils::storeImage("../data/AcmpK1.png", Image(AcmpK1), AcmpK1.rows(), AcmpK1.cols());
    Utils::storeImage("../data/AcmpK2.png", Image(AcmpK2), AcmpK2.rows(), AcmpK2.cols());

    /*
    TASK8
    Using Eigen create a black and white checkerboard image (as the one depicted below)
    with height and width equal to 200 pixels. Report the Euclidean norm of the matrix
    corresponding to the image.
    */
    int checkH = 200;
    int checkW = 200;
    int blockDim = 25;
    Image checkerboard = Image::Zero(checkH, checkW);
    Image whiteBlock = Image::Ones(blockDim, blockDim);

    int numBlocks = 8;
    int alternate = 0;
    for(int i = 0; i < numBlocks; ++i) {
        alternate = (i+1) % 2;
        for(int j = 0; j < numBlocks; j+=2) {
            checkerboard.block(i*blockDim,(j+alternate)*blockDim, blockDim, blockDim) = whiteBlock;
        }
    }

    // TEST
    Utils::storeImage("../data/checkerboard.png", checkerboard, checkerboard.rows(), checkerboard.cols());
    
    double checkerboardNorm = checkerboard.norm();
    std::cout << "[INFO] norm(checkerboard) = " << checkerboardNorm << std::endl;

    // TODO: check also without lambda expression
    /*
    TASK9
    Introduce a noise into the checkerboard image by adding random fluctuations of color
    ranging between [−50, 50] to each pixel. Export the resulting image in .png and upload it.
    */
    Image checkerboardNoisy = checkerboard.unaryExpr([](double val) -> double {
        double randValue = ((rand() % 101) - 50) / 255.0;

        if(val + randValue > 1.0) {
            return 1.0;
        }
        else if(val + randValue < 0.0) {
            return 0.0;
        }
        else {
            return val + randValue;
        }
    });

    Utils::storeImage("../data/checkerboardNoisy.png", checkerboardNoisy,
            checkerboardNoisy.rows(), checkerboardNoisy.cols());


    // FROM NOW ON: all the variables whose name starts with "checkerboard" are actually
    //  referring to the checkerboardNoisy image

    /*
    TASK10
    Using the SVD module of the Eigen library, perform a singular value decomposition of the
    matrix corresponding to the noisy image. Report the two largest computed singular values.
    */
    Eigen::BDCSVD<Image> svd2(checkerboardNoisy, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Vector checkerboardSingularValues = svd2.singularValues();
    Image checkerboardU = svd2.matrixU();
    Image checkerboardV = svd2.matrixV();
    Image checkerboardSigma = checkerboardSingularValues.asDiagonal();

    std::cout << "[INFO] σ1(checkerboardNoisy) = " << checkerboardSingularValues(0) << std::endl;
    std::cout << "[INFO] σ2(checkerboardNoisy) = " << checkerboardSingularValues(1) << std::endl;


    /*
    TASK11
    Starting from the previously computed SVD, create the matrices C and D defined in (1)
    assuming k = 5 and k = 10. Report the size of the matrices C and D.
    */
    double k3 = 5;
    Image Ck3 = checkerboardU.block(0,0,checkerboardU.rows(),k3);
    Image Dk3 = checkerboardV.block(0,0,checkerboardV.rows(),k3)*checkerboardSigma.block(0,0,k3,k3);

    double k4 = 10;
    Image Ck4 = checkerboardU.block(0,0,checkerboardU.rows(),k4);
    Image Dk4 = checkerboardV.block(0,0,checkerboardV.rows(),k4)*checkerboardSigma.block(0,0,k4,k4);

    std::cout << "[INFO] size(Ck3) (rows x cols) = " << Ck3.rows() << "x" << Ck3.cols() << std::endl;
    std::cout << "[INFO] size(Dk3) (rows x cols) = " << Dk3.rows() << "x" << Dk3.cols() << std::endl;
    std::cout << "[INFO] size(Ck4) (rows x cols) = " << Ck4.rows() << "x" << Ck4.cols() << std::endl;
    std::cout << "[INFO] size(Dk4) (rows x cols) = " << Dk4.rows() << "x" << Dk4.cols() << std::endl;

    /*
    TASK12
    Compute the compressed images as the matrix product CD^T (again for k = 5 and k = 10).
    Export and upload the resulting images in .png.
    */
    Image checkerbaordCmpK3 = Ck3*Dk3.transpose();
    Image checkerbaordCmpK4 = Ck4*Dk4.transpose();

    Utils::storeImage("../data/checkerbaordCmpK3.png", checkerbaordCmpK3,
            checkerbaordCmpK3.rows(), checkerbaordCmpK3.cols());
    Utils::storeImage("../data/checkerbaordCmpK4.png", checkerbaordCmpK4,
            checkerbaordCmpK4.rows(), checkerbaordCmpK4.cols());

    return 0;
}
