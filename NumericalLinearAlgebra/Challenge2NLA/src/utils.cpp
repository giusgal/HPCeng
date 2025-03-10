#include <iostream>
#include "../include/utils.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "../include/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../include/stb_image_write.h"

void Utils::loadImage(const char* path, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &image) {
    int w;
    int h;

    unsigned char* imageData = stbi_load(path, &w, &h, nullptr, 1);
    if (!imageData) {
        std::cerr << "[ERROR] couldn't load image" << std::endl;
        exit(1);
    }

    image.resize(h,w);
    for(int i = 0; i < h; ++i) {
        for(int j = 0; j < w; ++j) {
            image(i,j) = static_cast<double>(imageData[i*w+j]) / 255.0;
        }
    }
    stbi_image_free(imageData);
    std::cout << "[INFO] image " << path << " loaded" << std::endl;
}

void Utils::storeImage(const char* path, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const &image, int h, int w) {
    // casting + clipping
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp =
        image.unaryExpr([](double val) -> unsigned char {
            if(val > 1.0) {
                return 255;
            }
            if(val < 0.0) {
                return 0.0;
            }
            return static_cast<unsigned char>(val * 255.0);
        });

    if (stbi_write_png(path, w, h, 1, tmp.data(), w) == 0) {
        std::cerr << "[ERROR] couldn't store image" << std::endl;
        exit(1);
    }
    std::cout << "[INFO] image succesfully saved to " << path << std::endl;
}
