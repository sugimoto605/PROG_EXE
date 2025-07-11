#include <iostream>
#include <cmath>
#include <boost/timer/timer.hpp>                        //cpu-time
#include <boost/asio.hpp>                               //hostname
#include <boost/date_time/posix_time/posix_time.hpp>    //date-time
#include <boost/random.hpp>                             //random
#include <boost/math/special_functions/erf.hpp>         //ierf
#include <boost/multi_array.hpp>                        //multi-array

void boost_multi_array()
{
    boost::multi_array<double,3> A3(boost::extents[3][4][5]);
    std::cout << "A3 Dimension= " << A3.num_dimensions()
        << " Elements= " << A3.num_elements() << std::endl;
    for(int ndim=0;ndim<A3.num_dimensions();++ndim)
        std::cout << "A3 " << ndim << "-th max index=" << A3.shape()[ndim] << std::endl;
    A3[0][0][0] = 666.;
    A3[1][1][1] = 999.;
    A3[2][3][4] = -123;
    std::cout << "A3[0][0][0]=" << A3[0][0][0] << std::endl;
    std::cout << "A3[1][1][1]=" << A3[1][1][1] << std::endl;
    std::cout << "A3[2][3][4]=" << A3[2][3][4] << std::endl;
    for (int i = 0; i < A3.num_elements(); ++i)
        if (A3.data()[i]!=0) std::cout << "A3(" << i << ")=" << A3.data()[i] << std::endl;
    // index_bases()を追加すれば，添字の範囲を変更可能
    boost::array<boost::multi_array<double, 3>::index, 3> A3_base = {{5, 4, -5}};
    A3.reindex(A3_base);
    for (int ndim = 0; ndim < A3.num_dimensions(); ++ndim)
        std::cout << "A3 " << ndim << "-th max index=" << A3.shape()[ndim] 
        << " shifted to " << A3.index_bases()[ndim] << std::endl;
    std::cout << "A3[5][4][-5]=" << A3[5][4][-5] << std::endl;
    std::cout << "A3[6][5][-4]=" << A3[6][5][-4] << std::endl;
    std::cout << "A3[7][7][-1]=" << A3[7][7][-1] << std::endl;
    std::cout << "A3[8][8][0]=" << A3[8][8][0] << std::endl;
}

void boost_random()
{
    // 乱数生成器（メルセンヌ・ツイスタ）
    boost::random::mt19937 gen(static_cast<unsigned int>(std::time(0)));
    // 一様分布 [0,1)
    boost::random::uniform_real_distribution<> dist(0.0, 1.0);

    std::cout << "Boost random samples: ";
    for (int i = 0; i < 5; ++i) {
        std::cout << dist(gen) << " ";
    }
    std::cout << std::endl;
    // inverse-erf
    double x = 0.5;
    double y = boost::math::erf_inv(x); // erfの逆関数
    std::cout << "erf_inv(" << x << ") = " << y << std::endl;
}

void boost_timer()
{
    boost::timer::cpu_timer timer; // 計測開始

    // 計測したい処理
    double sum = 0;
    for (int i = 0; i < 100000000; ++i)
    {
        sum += std::sin(i * 0.000001);
    }
    boost::timer::cpu_times elapsed = timer.elapsed();
    std::cout << "user: " << elapsed.user / 1e9 << " sec" << std::endl;
    std::cout << "system: " << elapsed.system / 1e9 << " sec" << std::endl;
    std::cout << "wall: " << elapsed.wall / 1e9 << " sec" << std::endl;

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::asio::io_context io;
    std::string hostname = boost::asio::ip::host_name();
    std::cout << now << " on " << hostname << std::endl;
}

int main() {
    //boost_timer();
    //boost_random();
    boost_multi_array();
    return 0;
}