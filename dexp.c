// 近似誤差の計算が入ってない為、未完

#include "arb.h"

//関数f(x)=2*x^4
void func(arb_t res,const arb_t x,slong prec)
{
    arb_init(res);
    // res = x^4
    arb_pow_ui(res,x,4,prec);
    // res = 2*x^4
    arb_mul_ui(res,res,2,prec);
}

// double exponential fomula [a,b]
void dexp(arb_t res, slong a, slong b, slong n, slong prec)
{
    arb_init(res);
    // 刻み幅h
    arb_t h,ih,phi,phid,fx;
    arb_init(h);
    arb_init(ih);
    arb_init(phi);
    arb_init(phid);
    arb_init(fx);

    arb_t ba_sub,ba_add,pi,pi_sinh,cosh_ih;
    arb_init(ba_sub);
    arb_init(ba_add);
    arb_init(pi);
    arb_init(pi_sinh);
    arb_init(cosh_ih);

    //ba_sub = (b - a) * 0.5;
    arb_set_si(ba_sub,b);
    arb_sub_si(ba_sub,ba_sub,a,prec);
    arb_div_ui(ba_sub,ba_sub,2,prec);

    //ba_add = (b + a) * 0.5;
    arb_set_si(ba_add, b);
    arb_add_si(ba_add, ba_add, a, prec);
    arb_div_ui(ba_add, ba_add, 2, prec);

    // pi = π
    arb_const_pi(pi,prec);

    // h = log(3.0 * n) / n
    arb_set_si(h,n);
    arb_mul_ui(h,h,3,prec);
    arb_log(h,h,prec);
    arb_div_si(h,h,n,prec);



    for(int i = -n; i <= n; i++) {
        // ih = i * h
        arb_set(ih,h);
        arb_mul_si(ih,ih,i,prec);
        // pi_sinh = 0.5 * M_PI * sinh(i * h);
        arb_sinh(pi_sinh,ih,prec);
        arb_mul(pi_sinh,pi_sinh,pi,prec);
        arb_div_ui(pi_sinh,pi_sinh,2,prec);
        // phid = cosh(i*h)/(cosh(pi_sinh)^2)
        arb_cosh(phid,pi_sinh,prec);
        arb_inv(phid,phid,prec);
        arb_pow_ui(phid,phid,2,prec);
        arb_cosh(cosh_ih,ih,prec);
        arb_mul(phid,phid,cosh_ih,prec);
        // phi = ba_sub * tanh(pi_sinh) + ba_add
        arb_tanh(phi,pi_sinh,prec);
        arb_mul(phi,phi,ba_sub,prec);
        arb_add(phi,phi,ba_add,prec);
        // fx = func(phi)
        func(fx,phi,prec);
        // sum += fx*phid
        arb_mul(fx,fx,phid,prec);
        arb_add(res,res,fx,prec);
    }
    // sum *= 0.5 * M_PI * h * ba_sub
    arb_div_ui(res,res,2,prec);
    arb_mul(res,res,pi,prec);
    arb_mul(res, res, h, prec);
    arb_mul(res, res, ba_sub, prec);

    arb_clear(h);
    arb_clear(ih);
    arb_clear(phi);
    arb_clear(phid);
    arb_clear(fx);
    arb_clear(ba_sub);
    arb_clear(ba_add);
    arb_clear(pi);
    arb_clear(pi_sinh);
    arb_clear(cosh_ih);
}




int main()
{
    arb_t res, x;
    arb_init(res);
    arb_init(x);

    // 精度,表示桁数
    slong prec = 53, digit = 17;
    // 積分区間[a,b]
    slong a, b;
    a = -1;
    b = 1;
    // 分点数
    ulong n = 2;


    for(int i = 2; i<=64; i++) {
        dexp(res,a,b,n,prec);
        printf("[%5d]", i);
        arb_printn(res, digit, 0);
        flint_printf("\n");
        n++;
    }

    arb_clear(res);
    arb_clear(x);

    return 0;
}