#include <stdio.h>

const unsigned short ctab[] = {
#define Z(a) a, a + 1, a + 256, a + 257
#define Y(a) Z(a), Z(a + 2), Z(a + 512), Z(a + 514)
#define X(a) Y(a), Y(a + 4), Y(a + 1024), Y(a + 1028)
    X(0), X(8), X(2048), X(2056)
#undef X
#undef Y
#undef Z
};

const unsigned short utab[] = {
#define Z(a) 0x##a##0, 0x##a##1, 0x##a##4, 0x##a##5
#define Y(a) Z(a##0), Z(a##1), Z(a##4), Z(a##5)
#define X(a) Y(a##0), Y(a##1), Y(a##4), Y(a##5)
    X(0), X(1), X(4), X(5)
#undef X
#undef Y
#undef Z
};

int main(void)
{
    const int nitems = sizeof(utab) / sizeof(utab[0]);
    int i;

    for (i = 0; i < nitems; ++i)
    {
        printf("%5d, ", utab[i]);
        if (i > 0 && i % 5 == 0)
            puts("");
    }

    return 0;
}