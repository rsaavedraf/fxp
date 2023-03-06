/* SPDX-License-Identifier: MIT */
/*
 * bkm.c
 *
 * Initial implementations of lg2 and pow2 using long doubles.
 * Written to test the implemented algorithms in general,
 * before tailoring them for fxp's
 *
 * Also generates the tables of pre-calculated values for
 * long and int fxp implementations of BKM-based lg2 and pow2
 *
 * By Raul Saavedra, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "fxp_extern.h"

const long double FXP_ZERO_LD = 1.0E-124;

#define DASHES "=========================================\n"

static const long double A_2[] = {
    1.0000000000000000000000000000000000000000000000000000000000000000000000000000,
    0.5849625007211561814537389439478165087598144076924810604557526545410982276485,
    0.3219280948873623478703194294893901758648313930245806120547563958159347765589,
    0.1699250014423123629074778878956330175196288153849621209115053090821964552970,
    0.0874628412503394082540660108104043540112672823448206881266090643866965081686,
    0.0443941193584534376531019906736094674630459333742491317685543002674288465967,
    0.0223678130284545082671320837460849094932677948156179815932199216587899627785,
    0.0112272554232541203378805844158839407281095943600297940811823651462712311786,
    0.0056245491938781069198591026740666017211096815383520359072957784732489771013,
    0.0028150156070540381547362547502839489729507927389771959487826944878598909400,
    0.0014081943928083889066101665016890524233311715793462235597709051792834906001,
    0.0007042690112466432585379340422201964456668872087249334581924550139514213168,
    0.0003521774803010272377989609925281744988670304302127133979341729842842377649,
    0.0001760994864425060348637509459678580940163670081839283659942864068257522373,
    0.0000880524301221769086378699983597183301490534085738474534831071719854721939,
    0.0000440268868273167176441087067175806394819146645511899503059774914593663365,
    0.0000220136113603404964890728830697555571275493801909791504158295359319433723,
    0.0000110068476674814423006223021573490183469930819844945565597452748333526464,
    0.0000055034343306486037230640321058826431606183125807276574241540303833251704,
    0.0000027517197895612831123023958331509538486493412831626219340570294203116559,
    0.0000013758605508411382010566802834037147561973553922354232704569052932922954,
    0.0000006879304394358496786728937442939160483304056131990916985043387874690617,
    0.0000003439652607217645360118314743718005315334062644619363447395987584138324,
    0.0000001719826406118446361936972479533123619972434705828085978955697643547921,
    0.0000000859913228686632156462565208266682841603921494181830811515318381744650,
    0.0000000429956620750168703982940244684787907148132725669106053076409624949917,
    0.0000000214978311976797556164155504126645192380395989504741781512309853438587,
    0.0000000107489156388827085092095702361647949603617203979413516082280717515504,
    0.0000000053744578294520620044408178949217773318785601260677517784797554422804,
    0.0000000026872289172287079490026152352638891824761667284401180026908031182361,
    0.0000000013436144592400232123622589569799954658536700992739887706412976115422,
    0.0000000006718072297764289157920422846078078155859484240808550018085324187007,
    0.0000000003359036149273187853169587152657145221968468364663464125722491530858,
    0.0000000001679518074734354745159899223037458278711244127245990591908996412262,
    0.0000000000839759037391617577226571237484864917411614198675604731728132152582,
    0.0000000000419879518701918839775296677020135040214077417929807824842667285938,
    0.0000000000209939759352486932678195559552767641474249812845414125580747434389,
    0.0000000000104969879676625344536740142096218372850561859495065136990936290929,
    0.0000000000052484939838408141817781356260462777942148580518406975851213868092,
    0.0000000000026242469919227938296243586262369156865545638305682553644113887909,
    0.0000000000013121234959619935994960031017850191710121890821178731821983105443,
    0.0000000000006560617479811459709189576337295395590603644549624717910616347038,
    0.0000000000003280308739906102782522178545328259781415615142931952662153623493,
    0.0000000000001640154369953144623242936888032768768777422997704541618141646683,
    0.0000000000000820077184976595619616930350508356401599552034612281802599177300,
    0.0000000000000410038592488303636807330652208397742314215159774270270147020117,
    0.0000000000000205019296244153275153381695384157073687186580546938331088730952,
    0.0000000000000102509648122077001764119940017243502120046885379813510430378661,
    0.0000000000000051254824061038591928917243090559919209628584150482483994782302,
    0.0000000000000025627412030519318726172939815845367496027046030028595094737777,
    0.0000000000000012813706015259665053515049475574143952543145124550608158430592,
    0.0000000000000006406853007629833949364669629701200556369782295210193569318434,
    0.0000000000000003203426503814917330334121037829290364330169106716787999052925,
    0.0000000000000001601713251907458754080007074659337446341494733882570243497196,
    0.0000000000000000800856625953729399268240176265844257044861248416330071223615,
    0.0000000000000000400428312976864705191179247866966320469710511619971334577509,
    0.0000000000000000200214156488432353984854413866994246781519154793320684126179,
    0.0000000000000000100107078244216177339743404416874899847406043033792202127070,
    0.0000000000000000050053539122108088756700751579281894640362199287591340285355,
    0.0000000000000000025026769561054044400057638132352058574658089256646014899499,
    0.0000000000000000012513384780527022205455634651853807110362316427807660551208,
    0.0000000000000000006256692390263511104084521222346348012116229213309001913762,
    0.0000000000000000003128346195131755552381436585278035120438976487697544916191,
    0.0000000000000000001564173097565877776275512286165232838833090480508502328437,
    0.0000000000000000000782086548782938888158954641464170239072244145219054734086,
    0.0000000000000000000391043274391469444084776945327473574450334092075712154016,
    0.0000000000000000000195521637195734722043713378812583900953755962557525252782,
    0.0000000000000000000097760818597867361022187915943503728909029699365320287407,
    0.0000000000000000000048880409298933680511176764606054809062553340323879609794,
    0.0000000000000000000024440204649466840255609083961603140683286362962192177597,
    0.0000000000000000000012220102324733420127809717395445504379645613448652614939,
    0.0000000000000000000006110051162366710063906152551383735699323415812152114058,
    0.0000000000000000000003055025581183355031953399739107113727036860315024588989,
    0.0000000000000000000001527512790591677515976780735407368332862218276873443537,
    0.0000000000000000000000763756395295838757988410584167137033767056170417508383,
    0.0000000000000000000000381878197647919378994210346199431733717514843471513618,
    0.0000000000000000000000190939098823959689497106436628681671067254111334889005,
    0.0000000000000000000000095469549411979844748553534196582286585751228071408728,
    0.0000000000000000000000047734774705989922374276846068851506055906657137209047,
    0.0000000000000000000000023867387352994961187138442777065843718711089344045782,
    0.0000000000000000000000011933693676497480593569226324192944532044984865894525,
    0.0000000000000000000000005966846838248740296784614396011477934194852481410926,
    0.0000000000000000000000002983423419124370148392307506484490384140516252814304,
    0.0000000000000000000000001491711709562185074196153830361933046331030629430117,
    0.0000000000000000000000000745855854781092537098076934460888486730708440475045,
    0.0000000000000000000000000372927927390546268549038472050424734256652501673274,
    0.0000000000000000000000000186463963695273134274519237230207489851150821191330,
    0.0000000000000000000000000093231981847636567137259618916352525606281553180093,
    0.0000000000000000000000000046615990923818283568629809533488457973317312233323,
    0.0000000000000000000000000023307995461909141784314904785572277779202790023236,
    0.0000000000000000000000000011653997730954570892157452397493151087737428485431,
    0.0000000000000000000000000005826998865477285446078726199923328593402722606924,
    0.0000000000000000000000000002913499432738642723039363100255852559084863397344,
    0.0000000000000000000000000001456749716369321361519681550201473345138307215067,
    0.0000000000000000000000000000728374858184660680759840775119123438968122488047,
    0.0000000000000000000000000000364187429092330340379920387564158411083803465567,
    0.0000000000000000000000000000182093714546165170189960193783228378441837282509,
    0.0000000000000000000000000000091046857273082585094980096891901482445902524441,
    0.0000000000000000000000000000045523428636541292547490048446022564529197237262,
    0.0000000000000000000000000000022761714318270646273745024223029238091160103901
};

/*
 * Calculate the log2 of a long double number
 * (Implementing and testing the logarithm algorithm first
 * with floating points, before implementing it for fxp's)
 *
 * Based on the general algorithm to calculate binary logarithms
 * explained by Clay Turner in IEEE Signal Processing Magazine,
 * Sep/2010. His short article has D. E. Knuth's "The Art of
 * Computer Programming Vol 2: Seminumerical Algorithms",
 * 2nd ed., 1981 (pages 441 - 446) as the only reference.
 * Additional useful resources: CORDIC and/or BKM algorithms
 */
long double my_log2(long double x)
{
        if (x < 0.0) return FXP_UNDEF_LD;
        if (x <= FXP_ZERO_LD) return FXP_NINF_LD;
        long double z = x;
        int c = 0; // characteristic
        while (z >= 2.0) {
                c++;
                z /= 2;
        }
        while (z < 1.0) {
                c--;
                z *= 2;
        }
        // Here we have already calculated the log characteristic c, and
        // we have z satisfying: 1 <= z < 2, so we can use it to calculate
        // the mantissa
        //printf("log c: %d\n", c);
        long double m = 0.0; // mantissa
        long double b = 0.5;
        int nb = 64; // desired number of mantissa bits to process
        while (nb > 0) {
                z = z * z;
                if (z >= 2.0) {
                        z /= 2;
                        m += b;
                }
                b = b/2;
                nb--;
        }
        // Here we have already calculated the mantissa down to
        // nb bits of precision
        //printf("log m: %Lf\n", m);
        long double full_log = ((long double) c) + m;
        //printf("log(%Lf) = c + m = %Lf\n", x, full_log);
        return full_log;
}

/*
 * Implementation of the BKM L-Mode algorithm for log2 calculation
 * using long doubles:
 * https://en.wikipedia.org/wiki/BKM_algorithm
 */
long double my_log2_bkm(long double x, int nbits)
{
        if (x < 0.0) return FXP_UNDEF_LD;
        if (x <= FXP_ZERO_LD) return FXP_NINF_LD;
        long double z = x;
        int c = 0; // characteristic
        while (z >= 2.0) {
                c++;
                z /= 2;
        }
        while (z < 1.0) {
                c--;
                z *= 2;
        }
        // Here we have already calculated the log characteristic c, and
        // we have z satisfying: 1 <= z < 2
        long double xx = 1.0;
        long double yy = 0.0;
        long double ss = 0.5;
        // Notice ss starting with 0.5. Skipping the 1.0 because we know
        // zz < 2, so that first test if (zz < z) when ss == 1.0 would be
        // false for sure
        long double zz;
        for (int k = 1; k < nbits; k++) {
                // Notice that z is for sure in [1, 2), and yy
                // will remain in the range of log(z), and that is:
                // [0, 1), but zz can in fact get bigger than 2
                // repeatedly in these iterations, even if not by much.
                zz = xx + xx * ss;
                printf("\tlg2 bkm iteration %d, zz:%0.5LE\n", k, zz);
                if (zz <= z) {
                        xx = zz;
                        yy += A_2[k];
                        printf("\t\tUpdating yy:%0.5LE\n", yy);
                }
                ss *= 0.5;
        }
        long double full_log = ((long double) c) + yy;
        printf("log(%.5LE) = c (%d)+ m (%.5LE) = %.5LE\n", x, c, yy, full_log);
        return full_log;
}

/*
 * Implementation of the BKM E-Mode algorithm for 2^n calculation
 * using long doubles
 */
long double my_pow2_bkm(long double n, int nbits)
{
        // Notice that 2^n == 2^(whole(n) + frac(n)) ==
        // 2^whole(n) * 2^frac(n) = pow2(w) * pow2(f)
        long long w = truncl(n);
        long double frac = n - w;
        // Calculate pow2(w), and prepare argument for BKM
        long double pow2w;
        long double argument;
        if (n >= 0) {
                pow2w = (long double) (1l << w);
                argument = frac;
        } else {
                pow2w = 1.0 / ((long double) (1l << (-w + 1)));
                argument = 1 + frac; // Notice Argument is >= 0
        }
        printf("n:%.5Lf,  w:%lld,  argument:%.5LE\n", n, w, argument);
        // The BKM algorithm in E-Mode (for exponential)
        // calculates pow2(argument), argument in [0, 1) here
        long double x = 1.0, y = 0.0, s = 1.0;
        for (int k = 0; k < nbits; k++) {
                long double const  z = y + A_2[k];
                printf("\tpow2 bkm iteration %d, zz:%.5LE, x:%.5LE\n", \
                            k, z, x);
                if (z <= argument) {
                        y = z;
                        x = x + x*s;
                        printf("\t\tUpdating y:%.5LE x:%.5LE\n", y, x);
                }
                s *= 0.5;
        }
        // Here x == 2^argument == pow2(f)
        // Calculate the full 2^n = pow2w * pow2f
        long double full_pow2 = pow2w * x;
        printf("pow2(%.5Lf) = %.5Lf\n\n", n, full_pow2);
        return full_pow2;
}

int main(void)
{
        printf("\n%sGenerate table of log values for BKM.c\n%s", DASHES, DASHES);
        const int nfbits = 63;
        printf("\nValues for the BKM array, long hex format, %d frac bits:\n", nfbits);
        int n = sizeof(A_2) / sizeof(A_2[0]);
        int j = 1;
        unsigned long long mult = 1L << nfbits;
        printf("\nstatic const unsigned long FXP_BKM_LOGS_L[] = {\n");
        for (int i=0; i < n; i++) {
                long double numld = A_2[i];
                unsigned long long num = truncl( numld * mult );
                //printf("%d:\t0x%llX \t(%LE)\n", i, num, numld);
                printf("0x%llX", num);
                if (num == 0L) break;
                printf(", ");
                if (((i + 1) % 4) == 0) printf("\n");
        }
        printf("\n};\n");

        unsigned long long msbits, lsbits;
        printf("\nSplit values for the BKM arrays using only ints:\n");
        printf("\nFirst array:\n");
        printf("\nstatic const unsigned int FXP_BKM_LOGS[] = {\n");
        for (int i=0; i < n; i++) {
                long double numld = A_2[i];
                unsigned long long num = truncl( numld * mult );
                msbits = (num & 0xFFFFFFFF00000000) >> 32;
                printf("0x%llX", msbits);
                if (num == 0L) break;
                printf(", ");
                if (((i + 1) % 4) == 0) printf("\n");
        }
        printf("\n};\n");
        printf("\nXtra array:\n");
        printf("\nstatic const unsigned int FXP_BKM_LOGS_XTRA[] = {\n");
        for (int i=0; i < n; i++) {
                long double numld = A_2[i];
                unsigned long long num = truncl( numld * mult );
                lsbits = (num & 0x00000000FFFFFFFF);
                printf("0x%llX", lsbits);
                if (num == 0L) break;
                printf(", ");
                if (((i + 1) % 4) == 0) printf("\n");
        }
        printf("\n};\n");

        //long double x1 = my_log2_bkm(1.9999999999, 31);
        //long double x2 = my_log2_bkm(0.9999999999, 31);

        long double x3 = my_pow2_bkm(0, 31);
        //long double x4 = my_pow2_bkm(1.5, 31);
        //long double x5 = my_pow2_bkm(2.0, 31);
        //long double x6 = my_pow2_bkm(2.5, 31);

/*
        printf("\nPrecision of long double pow2() calculations using BKM:\n");
        long double invln2 = 1 / logl(2);
        for (int bkmd = 29; bkmd <= 31; bkmd++) {
                printf("\nChecking BKM depth %d\n", bkmd);
                // Making sure it works as expected even beyond 32 bits
                long double nums[] = {0.0, 0.25, 0.5, 1.0, 2.0, 3.1, 31.1, 40.1, \
                                    -0.25, -0.5, -1.0, -2.0, -3.1, -31.1, -40.1};
                int n = sizeof(nums) / sizeof(nums[0]);
                for (int i=0; i < n; i++) {
                        long double tgt = powl(2.0, nums[i]);
                        long double mypow2 = my_pow2_bkm(nums[i], bkmd);
                        printf("\n\tref: pow(2, %.3LE) = %.20LE\n",
                                nums[i], tgt);
                        printf("\tmy_pow2(%.3LE)     = %.20LE (- ref: %1.2LE)\n",
                                nums[i], mypow2, (mypow2 - tgt));
                }
        }
*/

}