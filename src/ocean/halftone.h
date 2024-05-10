#ifndef HALFTONE_H
#define HALFTONE_H


#include <vector>
class halftone
{
public:
    halftone();

private:
    int N = 0;
    int M = 0;
    float m_sigma = 24.f;
    std::vector<std::vector<int>> m_halftone;

};

#endif // HALFTONE_H
