#include "bashfest.h"

t_float lpp_ellipse(t_float x, LSTRUCT *eel, int nsects, t_float xnorm)
{
  register int m;
  t_float op;

  for(m=0;m<nsects;m++) {
    op = x + eel[m].c0 * eel[m].ps0 + eel[m].c2 * eel[m].ps1
      - eel[m].c1 * eel[m].ps2 - eel[m].c3 * eel[m].ps3;
    eel[m].ps1 = eel[m].ps0;
    eel[m].ps0 = x;
    eel[m].ps3 = eel[m].ps2;
    eel[m].ps2 = op;
    x = op;
  }
  return(x*xnorm);
}

void lpp_ellipset(t_float *list, LSTRUCT *eel, int  *nsects, t_float *xnorm)
{
/* the first argument in the list is the number of sections */
  int m,i;
  *nsects = (int)list[0];
  if(*nsects > MAXSECTS) {
    pd_error(0, "sorry, only configured for %d sections",MAXSECTS);
    return;
  }
  i=1;
  for(m=0;m<*nsects;m++) {
    eel[m].c0 = list[i++];
    eel[m].c1 = list[i++];
    eel[m].c2 = list[i++];
    eel[m].c3 = list[i++];
    eel[m].ps0 = eel[m].ps1 = eel[m].ps2 = eel[m].ps3 = 0;
  }
  *xnorm = list[i];
}
/*set biquad coefficients one time*/
void lpp_init_ellipse_data(t_float **a)
{
  /* 0: hipass at 200 */
  a[0][0] = 4;
  a[0][1] = 1.5156562;
  a[0][2] = -1.9958239;
  a[0][3] = 1;
  a[0][4] = 0.9965234;
  a[0][5] = -1.9996996;
  a[0][6] = 0.97229244;
  a[0][7] = 1;
  a[0][8] = 0.89313463;
  a[0][9] = 1.8828678;
  a[0][10] = -1.9561138;
  a[0][11] = 1;
  a[0][12] = 0.95845959;
  a[0][13] = -1.9999342;
  a[0][14] = -0.065893592;
  a[0][15] = 1;
  a[0][16] = 0.39565826;
  a[0][17] = 0.26225143;
  /* 1: hipass at 500 */
  a[1][0] = 3;
  a[1][1] = -1.9950633;
  a[1][2] = -1.9910746;
  a[1][3] = 1;
  a[1][4] = 0.99729187;
  a[1][5] = -1.9964182;
  a[1][6] = -1.9742933;
  a[1][7] = 1;
  a[1][8] = 0.98250468;
  a[1][9] = -1.999304;
  a[1][10] = -1.8149534;
  a[1][11] = 1;
  a[1][12] = 0.84353385;
  a[1][13] = 0.90419364;
  /* 2: bandpass 280 - 700 */
  a[2][0] = 4;
  a[2][1] = -1.9989934;
  a[2][2] = -1.9946771;
  a[2][3] = 1;
  a[2][4] = 0.99626146;
  a[2][5] = -1.9843098;
  a[2][6] = -1.9807532;
  a[2][7] = 1;
  a[2][8] = 0.99066977;
  a[2][9] = -1.9996779;
  a[2][10] = -1.9795816;
  a[2][11] = 1;
  a[2][12] = 0.9820447;
  a[2][13] = -1.9513627;
  a[2][14] = -1.965153;
  a[2][15] = 1;
  a[2][16] = 0.97142923;
  a[2][17] = 0.013949928;
  /* 3: lopass at 500 */
  a[3][0] = 3;
  a[3][1] = -1.9922858;
  a[3][2] = -1.9903447;
  a[3][3] = 1;
  a[3][4] = 0.99525722;
  a[3][5] = -1.9849712;
  a[3][6] = -1.9765264;
  a[3][7] = 1;
  a[3][8] = 0.97923558;
  a[3][9] = 1;
  a[3][10] = -0.98180316;
  a[3][11] = 0;
  a[3][12] = 0;
  a[3][13] = 0.0014021298;
  /* 4: lopass at 2K */
  a[4][0] = 3;
  a[4][1] = -1.9170388;
  a[4][2] = -1.9264647;
  a[4][3] = 1;
  a[4][4] = 0.99064223;
  a[4][5] = -1.8850187;
  a[4][6] = -1.9092573;
  a[4][7] = 1;
  a[4][8] = 0.95627234;
  a[4][9] = -1.4613313;
  a[4][10] = -1.8821271;
  a[4][11] = 0.99999996;
  a[4][12] = 0.89422169;
  a[4][13] = 0.0071020706;
  /* 5:
     f1,f2,f3= 400.0     900.0     2500.     ripple= 1.000     db= 40.00    */
  a[5][0] = 4;
  a[5][1] = -1.9594111;
  a[5][2] = -1.9900446;
  a[5][3] = 1;
  a[5][4] = 0.9932218;
  a[5][5] = -1.9986938;
  a[5][6] = -1.967975;
  a[5][7] = 1;
  a[5][8] = 0.98456911;
  a[5][9] = -1.8436241;
  a[5][10] = -1.9696535;
  a[5][11] = 1;
  a[5][12] = 0.9745592;
  a[5][13] = -1.9996708;
  a[5][14] = -1.9523041;
  a[5][15] = 1;
  a[5][16] = 0.96284515;
  a[5][17] = 0.0027629927;
  /* 6:  hipass at 500 */
  a[6][0] = 3;
  a[6][1] = -1.9950633;
  a[6][2] = -1.9910746;
  a[6][3] = 1;
  a[6][4] = 0.99729187;
  a[6][5] = -1.9964182;
  a[6][6] = -1.9742933;
  a[6][7] = 1;
  a[6][8] = 0.98250468;
  a[6][9] = -1.999304;
  a[6][10] = -1.8149534;
  a[6][11] = 1;
  a[6][12] = 0.84353385;
  a[6][13] = 0.90419364;
  /* 7: bandpass 1000-4000-6000 */
  a[7][0] = 4;
  a[7][1] = -1.9924766;
  a[7][2] = -1.9573893;
  a[7][3] = 1;
  a[7][4] = 0.97714717;
  a[7][5] = -1.2471186;
  a[7][6] = -1.6082873;
  a[7][7] = 1;
  a[7][8] = 0.91482231;
  a[7][9] = -1.9983103;
  a[7][10] = -1.8512918;
  a[7][11] = 1;
  a[7][12] = 0.88910013;
  a[7][13] = 0.033335148;
  a[7][14] = -1.6413378;
  a[7][15] = 0.99999998;
  a[7][16] = 0.78964041;
  a[7][17] = 0.0087452226;
  /* 8: hipass 4000-10000 */
  a[8][0] = 1;
  a[8][1] = -1.9896868;
  a[8][2] = -1.3953066;
  a[8][3] = 1;
  a[8][4] = 0.58943112;
  a[8][5] = 0.74811328;
  /* 9: bandpass 500-2500-3500 */
  a[9][0] = 6;
  a[9][1] = -1.9975736;
  a[9][2] = -1.9902167;
  a[9][3] = 1;
  a[9][4] = 0.99529287;
  a[9][5] = -1.7460823;
  a[9][6] = -1.853476;
  a[9][7] = 1;
  a[9][8] = 0.97721553;
  a[9][9] = -1.9984481;
  a[9][10] = -1.9737545;
  a[9][11] = 1;
  a[9][12] = 0.98056598;
  a[9][13] = -1.6166383;
  a[9][14] = -1.8408836;
  a[9][15] = 0.99999999;
  a[9][16] = 0.93097271;
  a[9][17] = -1.9997426;
  a[9][18] = -1.9320458;
  a[9][19] = 1;
  a[9][20] = 0.94629262;
  a[9][21] = -0.44018748;
  a[9][22] = -1.8664352;
  a[9][23] = 0.99999993;
  a[9][24] = 0.90871633;
  a[9][25] = 0.00044746789;
  /* 10: bp-300-400-2500-70dB */
  a[10][0] = 8;
  a[10][1] = -1.7823256;
  a[10][2] = -1.9938863;
  a[10][3] = 1;
  a[10][4] = 0.99712611;
  a[10][5] = -1.9981713;
  a[10][6] = -1.8579881;
  a[10][7] = 1;
  a[10][8] = 0.9825214;
  a[10][9] = -1.7151492;
  a[10][10] = -1.9844167;
  a[10][11] = 1;
  a[10][12] = 0.9884184;
  a[10][13] = -1.9986272;
  a[10][14] = -1.8447412;
  a[10][15] = 1;
  a[10][16] = 0.94374559;
  a[10][17] = -1.3382862;
  a[10][18] = -1.9602273;
  a[10][19] = 1;
  a[10][20] = 0.96717992;
  a[10][21] = -1.9994689;
  a[10][22] = -1.8529558;
  a[10][23] = 1;
  a[10][24] = 0.90889168;
  a[10][25] = 1;
  a[10][26] = -1.903171;
  a[10][27] = 0;
  a[10][28] = 0.92280038;
  a[10][29] = -1;
  a[10][30] = 0;
  a[10][31] = 0;
  a[10][32] = 0;
  a[10][33] = 0.00022546378;
}
