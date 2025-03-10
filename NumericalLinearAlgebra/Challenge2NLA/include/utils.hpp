#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Utils {

void loadImage(const char* path, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &image);
void storeImage(const char* path, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const &image, int h, int w);

};

#endif
