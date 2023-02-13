// 近似誤差の計算が入ってない、未完

#include "arb.h"

//関数f(x)=2*x^4
void func(arb_t res,const arb_t x,slong prec) {
    arb_init(res);
    // res = x^4
    arb_pow_ui(res,x,4,prec);
    // res = 2*x^4
    arb_mul_ui(res,res,2,prec);
}

// 複合台形公式 積分区間[a,b] 分点数n 精度prec
void trapezoidal(arb_t res,slong a,slong b,ulong n,slong prec)
{
    arb_init(res);
    // 刻み幅h
    arb_t h,x,fx,fa,fb,sum;
    arb_init(h); arb_init(x); arb_init(fx); arb_init(fa); arb_init(fb); arb_init(sum);

    // h = (b-a)/n
    arb_set_si(h,b);
    arb_sub_si(h,h,a,prec);
    arb_div_ui(h,h,n,prec);



    for(int i = 1; i < n; i++) {
        // x = a + h * i
        arb_mul_ui(x,h,i,prec);
        arb_add_si(x,x,a,prec);
        // fx = func(x)
        func(fx,x,prec);
        // sum += func(x)
        arb_add(sum,sum,fx,prec);
    }
    // arb_init(x);
    // fa = func(a)
    arb_set_si(x, a);
    func(fa,x,prec);
    // fb = func(b)
    arb_set_si(x, b);
    func(fb, x, prec);
    // res += 0.5 * (func(a) + func(b))
    arb_add(fx,fa,fb,prec);
    arb_div_ui(fx,fx,2,prec);
    arb_add(sum,sum,fx,prec);
    // res = h * sum
    arb_mul(res,h,sum,prec);

    arb_clear(h);
    arb_clear(x);
    arb_clear(fx);
    arb_clear(fa);
    arb_clear(fb);
    arb_clear(sum);
}


int main(void) {
    arb_t res,x;
    arb_init(res); arb_init(x);

    // 精度,表示桁数
    slong prec = 53 , digit = 17;
    // 積分区間[a,b]
    slong a,b;
    a = -1;
    b = 1;
    // 分点数
    ulong n = 1;

    for (int i = 1; i <= 65536; i *= 2) {
        trapezoidal(res, a, b, i, prec);
        printf("[%5d]",i);
        arb_printn(res, digit, 0);
        flint_printf("\n");
        n *= 2;
    }


    arb_clear(res); arb_clear(x);

    return 0;
}