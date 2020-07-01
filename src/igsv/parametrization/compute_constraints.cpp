#include <igsv/parametrization/compute_constraints.h>

#include <opencv2/core/eigen.hpp>
#include <opencv2/imgproc.hpp>

namespace IGSV {

  bool compute_constraints(const cv::Mat& _grey,               //
                           const cv::Mat& _bw,                 //
                           double _sw_avg,                     //
                           cv::Mat& _grey_blurred,             //
                           unsigned& _n_dark,                  //
                           Eigen::MatrixXd& _pix_weight,       //
                           Eigen::MatrixXi& _pix2band,         //
                           Eigen::VectorXi& _PX_I,             //
                           Eigen::VectorXi& _PX_J,             //
                           Eigen::VectorXcd& _PX_XY,           //
                           Eigen::VectorXcd& _PX_SobelTangent, //
                           Eigen::VectorXd& _PX_SobelWeight) {

    using namespace cv;

    // blur
    const double blur_sigma = 0.1;
    cv::GaussianBlur(_grey, _grey_blurred, cv::Size(5, 5), blur_sigma);

    Mat dx, dy;
    Mat dx2, dy2, dxy;
    Sobel(_grey_blurred, dx, CV_64F, 1, 0);
    Sobel(_grey_blurred, dy, CV_64F, 0, 1);
    dy = -dy; //// !!! FLIP Y COORD !!!
    multiply(dx, dx, dx2);
    multiply(dy, dy, dy2);
    multiply(dx, dy, dxy);

    GaussianBlur(dx2, dx2, Size(5, 5), 1.0);
    GaussianBlur(dy2, dy2, Size(5, 5), 1.0);
    GaussianBlur(dxy, dxy, Size(5, 5), 1.0);

    // OpenCV --> Eigen
    Mat planes[] = { dx, dy };
    Mat cvG;
    merge(planes, 2, cvG);
    Eigen::MatrixXcd grad;
    cv2eigen(cvG, grad);

    Eigen::MatrixXcd tauTimesGmag = grad * std::complex<double>(0.0, 1.0);
    Eigen::MatrixXd gMag          = tauTimesGmag.cwiseAbs();

    double maxGradMag = gMag.maxCoeff();
    for (int i = 0; i < _bw.rows; ++i)
      for (int j = 0; j < _bw.cols; ++j)
        if (fabs(gMag(i, j) / maxGradMag) < 0.1) {
          gMag(i, j)         = 0;
          tauTimesGmag(i, j) = 0;
        }

    Eigen::MatrixXd gMagNoZeros = gMag;
    for (int i = 0; i < _bw.rows; ++i)
      for (int j = 0; j < _bw.cols; ++j)
        if (fabs(gMag(i, j)) < 1e-10)
          gMagNoZeros(i, j) = 1;

    // compute and store tau
    const Eigen::MatrixXcd tau           = tauTimesGmag.array() / gMagNoZeros.array();
    const Eigen::MatrixXcd tauTimesGmag2 = tauTimesGmag.array().pow(2);

    // calculate laplacian
    // TODO : get rid of so many copies back and forth
    Mat tauTimesGmag2Re, tauTimesGmag2Im;
    const Eigen::MatrixXd x = tauTimesGmag2.real();
    const Eigen::MatrixXd y = tauTimesGmag2.imag();
    eigen2cv(x, tauTimesGmag2Re);
    eigen2cv(y, tauTimesGmag2Im);

    Mat kernel              = Mat::ones(3, 3, CV_64F);
    kernel.at<double>(1, 1) = 0;

    Mat Lx, Ly;
    filter2D(tauTimesGmag2Re, Lx, -1, kernel);
    filter2D(tauTimesGmag2Im, Ly, -1, kernel);

    Eigen::MatrixXd Lx_eig, Ly_eig;
    cv2eigen(Lx, Lx_eig);
    cv2eigen(Ly, Ly_eig);

    Eigen::MatrixXcd mse = Lx_eig.cast<std::complex<double>>() //
                           + Ly_eig.cast<std::complex<double>>() * std::complex<double>(0, 1);
    Eigen::MatrixXd mseNorm = mse.cwiseAbs();
    for (int i = 0; i < _bw.rows; ++i)
      for (int j = 0; j < _bw.cols; ++j)
        if (mseNorm(i, j) < 1e-10)
          mseNorm(i, j) = 1;

    const Eigen::MatrixXcd tau2 = tau.array().pow(2);

    mse = mse.array() / mseNorm.array();
    mse = mse - tau2;

    for (int i = 0; i < _bw.rows; ++i)
      for (int j = 0; j < _bw.cols; ++j)
        if (fabs(gMag(i, j)) < 1e-10) {
          mse(i, j) = 0;
        }

    Eigen::MatrixXd weight = mse.cwiseAbs();
    weight                 = Eigen::MatrixXd::Ones(_bw.rows, _bw.cols) - weight / weight.maxCoeff();
    for (int i = 0; i < _bw.rows; ++i)
      for (int j = 0; j < _bw.cols; ++j)
        if (fabs(gMag(i, j)) < 1e-10)
          weight(i, j) = 0;

    // store weights for all pixels (used for smoothness term of frame field energy)
    _pix_weight = weight;

    // store narrow band data
    _n_dark = _bw.cols * _bw.rows - countNonZero(_bw);
    _pix2band.resize(_bw.rows, _bw.cols);
    _pix2band.fill(-1);
    _PX_I.resize(_n_dark, 1);
    _PX_J.resize(_n_dark, 1);
    _PX_XY.resize(_n_dark, 1);
    _PX_SobelTangent.resize(_n_dark, 1);
    _PX_SobelWeight.resize(_n_dark, 1);

    unsigned k = 0;
    for (int i = 0; i < _bw.rows; i++)
      for (int j = 0; j < _bw.cols; j++)
        if (_bw.at<unsigned char>(i, j) == 0) {
          _pix2band(i, j)     = k;
          _PX_I(k)            = i;
          _PX_J(k)            = j;
          _PX_XY(k)           = std::complex<double>((double)(j) + 0.5, (double)(_bw.rows - 1 - i) + 0.5);
          _PX_SobelTangent(k) = tau(i, j);
          _PX_SobelWeight(k)  = _pix_weight(i, j);
          k++;
        }

    // std::cout << "  " << _n_dark << " black pixels" << std::endl;

    return true;
  }

} // namespace IGSV
