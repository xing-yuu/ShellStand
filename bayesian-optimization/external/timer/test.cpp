#include <random>
#include <timer.hpp>

void perform_long_calculations()
{
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::uniform_real_distribution<> dist(-1.0, 1.0);
    
    for (int i = 0; i < 10000000; ++ i)
    {
        const double value = dist(engine) + dist(engine);
    }
}

int main()
{
    timer::Timer timer;
    perform_long_calculations();    
    return 0;
}
