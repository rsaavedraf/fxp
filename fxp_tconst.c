/*
 * fxp_tconst.c
 *
 * Auxiliary program to generate fxp versions of some
 * important/transcendental constants
 *
 * By Raul Saavedra, Bonn Germany
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include "fxp.h"
#include "fxp_aux.h"
#include "fxp_tconst.h"


#define DASHES "========================\n"

#define BKM_DEPTH 32

#define STR_TEST1 "0000000000000000011111100000000000"
#define STR_TEST2 "0000000000000000100100000000000000"
#define STR_TEST3 "0000000000000001111100000000000000"

// Euler's constant (e) in decimal and binary
#define STR_E_DEC  "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274"
#define STR_E_BIN "10.1011011111100001010100010110001010001010111011010010101001101010101111110111000101011000100000001001"

// Pi in decimal and binary
#define STR_PI_DEC  "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
#define STR_PI_BIN  "11.0010010000111111011010101000100010000101101000110000100011010011000100110001100110001010001011100000"

// ln(2)
#define STR_LN_2_DEC ".693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687"

// log10(2)
#define STR_LG10_2_DEC ".301029995663981195213738894724493026768189881462108541310427461127108189274424509486927252118186172040684"

// log2(e)
//#define STR_LG2_E_DEC   "1.4426950408889634073599246810018921374266459541529859341354494069311092191811850798855266228935063444"

// log2(10)
//#define STR_LG2_10_DEC "3.32192809488736234787031942948939017586483139302458061205475639581593477660862521585013974335937015"



// Sqrt(2)
//#define STR_SQRT2_BIN "1.01101010000010011110011001100111111100111011110011001001000010001011001011111011000100110110011011"
//#define STR_SQRT2_DEC "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157"

#define ZERO 1.0E-124

// Decimal values from table A_2 in:
// https://en.wikipedia.org/wiki/BKM_algorithm
// The values here are simply: a[k] = log2(1 + 0.5^k)
// Their sum can add up to 2.25+
static char * str_A_2[] = {
    // 76 decimal digits, 100 numbers in total
    //0----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*-
    "1.0000000000000000000000000000000000000000000000000000000000000000000000000000",
    "0.5849625007211561814537389439478165087598144076924810604557526545410982276485",
    "0.3219280948873623478703194294893901758648313930245806120547563958159347765589",
    "0.1699250014423123629074778878956330175196288153849621209115053090821964552970",
    "0.0874628412503394082540660108104043540112672823448206881266090643866965081686",
    "0.0443941193584534376531019906736094674630459333742491317685543002674288465967",
    "0.0223678130284545082671320837460849094932677948156179815932199216587899627785",
    "0.0112272554232541203378805844158839407281095943600297940811823651462712311786",
    "0.0056245491938781069198591026740666017211096815383520359072957784732489771013",
    "0.0028150156070540381547362547502839489729507927389771959487826944878598909400",
    "0.0014081943928083889066101665016890524233311715793462235597709051792834906001",
    "0.0007042690112466432585379340422201964456668872087249334581924550139514213168",
    "0.0003521774803010272377989609925281744988670304302127133979341729842842377649",
    "0.0001760994864425060348637509459678580940163670081839283659942864068257522373",
    "0.0000880524301221769086378699983597183301490534085738474534831071719854721939",
    "0.0000440268868273167176441087067175806394819146645511899503059774914593663365",
    "0.0000220136113603404964890728830697555571275493801909791504158295359319433723",
    "0.0000110068476674814423006223021573490183469930819844945565597452748333526464",
    "0.0000055034343306486037230640321058826431606183125807276574241540303833251704",
    "0.0000027517197895612831123023958331509538486493412831626219340570294203116559",
    "0.0000013758605508411382010566802834037147561973553922354232704569052932922954",
    "0.0000006879304394358496786728937442939160483304056131990916985043387874690617",
    "0.0000003439652607217645360118314743718005315334062644619363447395987584138324",
    "0.0000001719826406118446361936972479533123619972434705828085978955697643547921",
    "0.0000000859913228686632156462565208266682841603921494181830811515318381744650",
    "0.0000000429956620750168703982940244684787907148132725669106053076409624949917",
    "0.0000000214978311976797556164155504126645192380395989504741781512309853438587",
    "0.0000000107489156388827085092095702361647949603617203979413516082280717515504",
    "0.0000000053744578294520620044408178949217773318785601260677517784797554422804",
    "0.0000000026872289172287079490026152352638891824761667284401180026908031182361",
    "0.0000000013436144592400232123622589569799954658536700992739887706412976115422",
    "0.0000000006718072297764289157920422846078078155859484240808550018085324187007",
    "0.0000000003359036149273187853169587152657145221968468364663464125722491530858",
    "0.0000000001679518074734354745159899223037458278711244127245990591908996412262",
    "0.0000000000839759037391617577226571237484864917411614198675604731728132152582",
    "0.0000000000419879518701918839775296677020135040214077417929807824842667285938",
    "0.0000000000209939759352486932678195559552767641474249812845414125580747434389",
    "0.0000000000104969879676625344536740142096218372850561859495065136990936290929",
    "0.0000000000052484939838408141817781356260462777942148580518406975851213868092",
    "0.0000000000026242469919227938296243586262369156865545638305682553644113887909",
    "0.0000000000013121234959619935994960031017850191710121890821178731821983105443",
    "0.0000000000006560617479811459709189576337295395590603644549624717910616347038",
    "0.0000000000003280308739906102782522178545328259781415615142931952662153623493",
    "0.0000000000001640154369953144623242936888032768768777422997704541618141646683",
    "0.0000000000000820077184976595619616930350508356401599552034612281802599177300",
    "0.0000000000000410038592488303636807330652208397742314215159774270270147020117",
    "0.0000000000000205019296244153275153381695384157073687186580546938331088730952",
    "0.0000000000000102509648122077001764119940017243502120046885379813510430378661",
    "0.0000000000000051254824061038591928917243090559919209628584150482483994782302",
    "0.0000000000000025627412030519318726172939815845367496027046030028595094737777",
    "0.0000000000000012813706015259665053515049475574143952543145124550608158430592",
    "0.0000000000000006406853007629833949364669629701200556369782295210193569318434",
    "0.0000000000000003203426503814917330334121037829290364330169106716787999052925",
    "0.0000000000000001601713251907458754080007074659337446341494733882570243497196",
    "0.0000000000000000800856625953729399268240176265844257044861248416330071223615",
    "0.0000000000000000400428312976864705191179247866966320469710511619971334577509",
    "0.0000000000000000200214156488432353984854413866994246781519154793320684126179",
    "0.0000000000000000100107078244216177339743404416874899847406043033792202127070",
    "0.0000000000000000050053539122108088756700751579281894640362199287591340285355",
    "0.0000000000000000025026769561054044400057638132352058574658089256646014899499",
    "0.0000000000000000012513384780527022205455634651853807110362316427807660551208",
    "0.0000000000000000006256692390263511104084521222346348012116229213309001913762",
    "0.0000000000000000003128346195131755552381436585278035120438976487697544916191",
    "0.0000000000000000001564173097565877776275512286165232838833090480508502328437",
    "0.0000000000000000000782086548782938888158954641464170239072244145219054734086",
    "0.0000000000000000000391043274391469444084776945327473574450334092075712154016",
    "0.0000000000000000000195521637195734722043713378812583900953755962557525252782",
    "0.0000000000000000000097760818597867361022187915943503728909029699365320287407",
    "0.0000000000000000000048880409298933680511176764606054809062553340323879609794",
    "0.0000000000000000000024440204649466840255609083961603140683286362962192177597",
    "0.0000000000000000000012220102324733420127809717395445504379645613448652614939",
    "0.0000000000000000000006110051162366710063906152551383735699323415812152114058",
    "0.0000000000000000000003055025581183355031953399739107113727036860315024588989",
    "0.0000000000000000000001527512790591677515976780735407368332862218276873443537",
    "0.0000000000000000000000763756395295838757988410584167137033767056170417508383",
    "0.0000000000000000000000381878197647919378994210346199431733717514843471513618",
    "0.0000000000000000000000190939098823959689497106436628681671067254111334889005",
    "0.0000000000000000000000095469549411979844748553534196582286585751228071408728",
    "0.0000000000000000000000047734774705989922374276846068851506055906657137209047",
    "0.0000000000000000000000023867387352994961187138442777065843718711089344045782",
    "0.0000000000000000000000011933693676497480593569226324192944532044984865894525",
    "0.0000000000000000000000005966846838248740296784614396011477934194852481410926",
    "0.0000000000000000000000002983423419124370148392307506484490384140516252814304",
    "0.0000000000000000000000001491711709562185074196153830361933046331030629430117",
    "0.0000000000000000000000000745855854781092537098076934460888486730708440475045",
    "0.0000000000000000000000000372927927390546268549038472050424734256652501673274",
    "0.0000000000000000000000000186463963695273134274519237230207489851150821191330",
    "0.0000000000000000000000000093231981847636567137259618916352525606281553180093",
    "0.0000000000000000000000000046615990923818283568629809533488457973317312233323",
    "0.0000000000000000000000000023307995461909141784314904785572277779202790023236",
    "0.0000000000000000000000000011653997730954570892157452397493151087737428485431",
    "0.0000000000000000000000000005826998865477285446078726199923328593402722606924",
    "0.0000000000000000000000000002913499432738642723039363100255852559084863397344",
    "0.0000000000000000000000000001456749716369321361519681550201473345138307215067",
    "0.0000000000000000000000000000728374858184660680759840775119123438968122488047",
    "0.0000000000000000000000000000364187429092330340379920387564158411083803465567",
    "0.0000000000000000000000000000182093714546165170189960193783228378441837282509",
    "0.0000000000000000000000000000091046857273082585094980096891901482445902524441",
    "0.0000000000000000000000000000045523428636541292547490048446022564529197237262",
    "0.0000000000000000000000000000022761714318270646273745024223029238091160103901"
};

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

static const unsigned long BKM_HEX_LOGS[] = {
    0x4000000000000000, 0x2570068E7EF5A27D, 0x149A784BCD1B8B50, 0xAE00D1CFDEB43FB,
    0x598FDBEB244C5B5, 0x2D75A6EB1DFB0F1, 0x16E79685C2D229E, 0xB7F285B778428E,
    0x5C2711B5EAB1DE, 0x2E1F07FE14EACA, 0x1712653743F454, 0xB89EB17BCABE1,
    0x5C523B0A86FF2, 0x2E29D623F4A6C, 0x1715193B17D35, 0xB8A982801725,
    0x5C54EF6A3E08, 0x2E2A833FB72C, 0x171544828311, 0xB8AA2F9EB95,
    0x5C551AB2053, 0x2E2A8E11ACC, 0x1715473700F, 0xB8AA3A70B1,
    0x5C551D6683, 0x2E2A8EBECC, 0x1715476248, 0xB8AA3B1DD,
    0x5C551D91C, 0x2E2A8EC99, 0x17154764F, 0xB8AA3B28,
    0x5C551D94, 0x2E2A8ECA, 0x17154765, 0xB8AA3B2,
    0x5C551D9, 0x2E2A8EC, 0x1715476, 0xB8AA3B,      // <---- *
    0x5C551D, 0x2E2A8E, 0x171547, 0xB8AA3,
    0x5C551, 0x2E2A8, 0x17154, 0xB8AA,
    0x5C55, 0x2E2A, 0x1715, 0xB8A,
    0x5C5, 0x2E2, 0x171, 0xB8,
    0x5C, 0x2E, 0x17, 0xB,
    0x5, 0x2, 0x1
    // Starting with the row marked with the *, each entry is exactly
    // a 4 -bit right-shift of the value 4 positions earlier
};


/*
 * Return binary expansion (up to nbits) of a string representation
 * of a binary number, with last bit either verbatim or rounded
 */
unsigned long long bex_from_bin(char * pbinnum, int wbits, int fbits, int rounded)
{
        char * p = pbinnum;
        if ((wbits <= 0) || (wbits > 31)) wbits = 31;
        if ((fbits <= 0) || (fbits + wbits > 63)) fbits = 63 - wbits;
        int roomleft = wbits + fbits;
        unsigned long long bnum = 0;
        while ((roomleft > 0) && (*p != '\0')) {
            char cbit = *p;
            if (cbit != '.') {
                bnum = (bnum << 1) | ((cbit == '1')? 1: 0);
                roomleft--;
            }
            p++;
        }
        if (roomleft > 0) {
            bnum = bnum << roomleft;
            // Source string with too few bits to need any rounding,
            // so return already
            return bnum;
        }
        if (rounded) {
            //printf("  (Next bit:");
            int nextbit = 0;
            char cbit = *p;
            //printf("%c", cbit);
            int rbit = ((cbit  == '1')? 1: 0);
            //printf(") -> adding %d", rbit);
            bnum += rbit;
        }
        return bnum;
}

/*
 * Return binary expansion (up to wbits.fbits) of a string representation
 * of a decimal number, with last bit either verbatim or rounded
 */
unsigned long long bex_from_dec(char * pdecnum, int wbits, int fbits, int rounded)
{
        char * p = pdecnum;
        if ((wbits < 0) || (wbits > 31)) wbits = 31;
        if ((fbits < 0) || (fbits + wbits > 64)) fbits = 63 - wbits;
        // Process whole part first
        int nwd = 0;
        int ptens = 1;
        while (*p != '.') {
            p++;
            ptens *= 10;
            nwd++;
        }
        int wnum = 0;
        p = pdecnum;
        while (*p != '.') {
            char cdigit = *p;
            ptens /= 10;
            wnum += ((int) cdigit - '0') * ptens;
            p++;
        }
        //printf("Whole part is d:%d (%d dec digits)\n", wnum, nwd);
        long long binnum = 0;
        int nb = 0;
        while (wnum > 0) {
            int r = wnum % 2;
            //printf("Appending a %d to the whole part\n", r);
            binnum = (r << nb) | binnum;
            nb++; // here count the whole bits
            wnum /= 2;
        }
        //printf("binnum: (x:%llx. , %d bits)\n", binnum, nb);

        // Process fractional part
        p++; // skip the '.'
        int ndd = nwd;
        long double ftens = 1.0/10.0;
        long double frac = 0;
        int ndig = 0;
        while (*p != '\0') { //&& (ftens > ZERO)) {
            char cdigit = *p;
            frac += ((int) cdigit - '0') * ftens;
            ftens /= 10;
            ndig++;
            ndd++;
            p++;
        }
        //printf("From %d decimal frac digits, fraction is: %1.20Lf \n", ndig, frac);
        int roomleft = fbits;
        long double m;
        while (roomleft > 0) {
            m = frac * 2.0;
            if (m >= 1.0) {
                binnum = ((binnum << 1) | 1);
                frac = m - 1.0;
                //printf("1");
            } else {
                binnum = (binnum << 1);
                frac = m;
                //printf("0");
            }
            roomleft--;
        }
        if (rounded) {
            m = frac * 2.0;
            int rbit = ((m >= 1.0)? 1: 0);
            binnum += rbit;
            //printf(" (+ round bit %d)", rbit);
        }
        //printf("\n");
        return binnum;
}

/*
 * Return a long double from a str representation of a decimal number
 */
long double Lf_from_dec(char * pdecnum)
{
        char * p = pdecnum;
        // Process whole part first
        int nwd = 0;
        int ptens = 1;
        while (*p != '.') {
            p++;
            ptens *= 10;
            nwd++;
        }
        long double num = 0.0;
        p = pdecnum;
        while (*p != '.') {
            char cdigit = *p;
            ptens /= 10;
            num += ((int) cdigit - '0') * ptens;
            p++;
        }

        // Process fractional part
        p++; // skip the '.'
        long double tenth = 1.0 / 10.0;
        long double ftens = tenth;
        while (*p != '\0') {
            char cdigit = *p;
            num += ((int) cdigit - '0') * ftens;
            ftens *= tenth;
            p++;
        }
        return num;
}

long double Lf_from_bin(char * pbinnum)
{
        char * p = pbinnum;
        // Process whole part first
        int nwd = 0;
        int ptwo = 1;
        while (*p != '.') {
            p++;
            ptwo *= 2;
            nwd++;
        }
        long double num = 0.0;
        p = pbinnum;
        while (*p != '.') {
            char cdigit = *p;
            ptwo /= 2;
            num += ((int) cdigit - '0') * ptwo;
            p++;
        }

        // Process fractional part
        p++; // skip the '.'
        long double half = 0.5;
        long double frac = half;
        while (*p != '\0') {
            char cdigit = *p;
            num += ((int) cdigit - '0') * frac;
            frac *= half;
            p++;
        }
        return num;
}

#define LONGLONG_BITS ((int) (sizeof(long long) * 8))

// Calculate the log2 of a long double number
// (Implementing and testing the logarithm algorithm first
// with floating points, before implementing it for fxp's)
//
// Based on the general algorithm to calculate binary logarithms
// explained by Clay Turner in IEEE Signal Processing Magazine,
// Sep/2010. His short article has D. E. Knuth's "The Art of
// Computer Programming Vol 2: Seminumerical Algorithms",
// 2nd ed., 1981 (pages 441 - 446) as the only reference.
// Additional useful resources: CORDIC and/or BKM algorithms
long double my_log2(long double x)
{
        if (x < 0.0) return FXP_UNDEF_LD;
        if (x <= ZERO) return FXP_NINF_LD;
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

// Implementation of the bkm algorithm using long doubles:
// https://en.wikipedia.org/wiki/BKM_algorithm
long double my_log2_bkm(long double x, int nbits)
{
        if (x < 0.0) return FXP_UNDEF_LD;
        if (x <= ZERO) return FXP_NINF_LD;
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
        // Notice ss starting with 0.5 because our first table entry is not
        // the 1.00 but the next one. Skipping the 1.0 because we know
        // zz < 2, so that first test if (zz < z) when ss == 1.0 would be
        // false for sure
        long double zz;
        for (int k = 1; k < nbits; k++) {
                zz = xx + xx * ss;
                if (zz <= z) {
                        xx = zz;
                        yy += A_2[k];
                }
                ss *= 0.5;
        }
        long double full_log = ((long double) c) + yy;
        //printf("log(%Lf) = c + m = %Lf\n", x, full_log);
        return full_log;
}

int main(void)
{
        printf("\n%sfxp_tconst.c\n%s", DASHES, DASHES);
        print_sys_info();

        //printf("\ntest1 %llx\n", bex_from_bin(STR_TEST1, 16, 1));
        //printf("\ntest2 %llx\n", bex_from_bin(STR_TEST2, 16, 1));
        //printf("\ntest3 %llx\n", bex_from_bin(STR_TEST3, 16, 1));

        printf("\nPrecisions of e on this system:\n");
        printf("                     ----*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8----*----9----*----0\n");
        const long double E_EXP = (long double) exp(1);
        const long double MY_E_DEC = Lf_from_dec(STR_E_DEC);
        const long double MY_E_BIN = Lf_from_bin(STR_E_BIN);
        const long double E_EXPL = (long double) expl(1);
        printf("e from M_E       : %.64LE (- expl(): %1.4LE)\n",
                (long double) M_E, (M_E - E_EXPL));
        printf("e from exp()     : %.64LE (- expl(): %1.4LE)\n",
                E_EXP, (E_EXP - E_EXPL));
        printf("My e from dec str: %.64LE (- expl(): %1.4LE)\n",
                MY_E_DEC, (MY_E_DEC - E_EXPL));
        printf("My e from bin str: %.64LE (- expl(): %1.4LE)\n",
                MY_E_BIN, (MY_E_BIN - E_EXPL));
        printf("e from expl()    : %.64LE (- expl(): %1.4LE)\n",
                E_EXPL, (E_EXPL - E_EXPL));
        printf("'True' e         : %s\n\n", STR_E_DEC);

        printf("e as binary: %s\n", STR_E_BIN);
        printf("e as fxp (31 bits) %llx\n", bex_from_bin(STR_E_BIN, 2, 29, 1));
        printf("e as fxp (63 bits) %llx\n\n", bex_from_bin(STR_E_BIN, 2, 61, 1));

        printf("pi as binary: %s\n", STR_PI_BIN);
        printf("pi as fxp (31 bits) %llx\n", bex_from_bin(STR_PI_BIN, 2, 29, 1));
        printf("pi as fxp (63 bits) %llx\n\n", bex_from_bin(STR_PI_BIN, 2, 61, 1));

        printf("ln(2) as decimal: %s\n", STR_LN_2_DEC);
        printf("ln(2) as fxp (31 bits) %llx\n", bex_from_dec(STR_LN_2_DEC, 0, 31, 1));
        printf("ln(2) as fxp (63 bits) %llx\n\n", bex_from_dec(STR_LN_2_DEC, 0, 63, 1));

        printf("lg10(2) as decimal: %s\n", STR_LG10_2_DEC);
        printf("lg10(2) as fxp (31 bits) %llx\n", bex_from_dec(STR_LG10_2_DEC, 0, 31, 1));
        printf("lg10(2) as fxp (63 bits) %llx\n\n", bex_from_dec(STR_LG10_2_DEC, 0, 63, 1));

        // Generate fxp's for the logs table of numbers to be used
        // by BKM algorithm

        /*
        int nfbits = 62;
        printf("\nPrinting values for the BKM arrays, in long hex format:\n");
        fxp_set_frac_bits(nfbits);
        int n = sizeof(str_A_2) / sizeof(str_A_2[0]);
        int j = 1;
        long double sumlogs2 = 0;
        for (int i=0; i < n; i++) {
            unsigned long long num = bex_from_dec(str_A_2[i], 0, nfbits, 0);
            long double numld = Lf_from_dec(str_A_2[i]);
            //printf("%d:\t0x%llX \t(%LE)\n", i, num, numld);
            sumlogs2 += numld;
            // The full long:
            //printf("%16llX, \n", num);
            unsigned long long msb, lsb;
            msb = (num & 0xFFFFFFFF00000000) >> 32;
            lsb = (num & 0x00000000FFFFFFFF);
            // bytes in the long
            //printf("0x%llX, \n", msb);
            printf("0x%8llX, \n", lsb);
            //printf("%8llX-%8llX\n", msb, lsb);
            //if (j % 4 == 0) printf("\n");
            j++;
        }
        //printf("Sum of values in log table: %1.5LE\n", sumlogs2);
        */

        /*
        printf("\nPrecision of long double log2 calculations using BKM:\n");
        long double invln2 = 1 / logl(2);
        for (int bkmd = 27; bkmd <= 35; bkmd++) {
                printf("\nChecking BKM depth %d\n", bkmd);
                long double nums[] = {0.0001, 0.00000001, 0.00000000001};
                int n = sizeof(nums) / sizeof(nums[0]);
                for (int i=0; i < n; i++) {
                        long double tgt = logl(nums[i]) * invln2;
                        long double mylg2 = my_log2(nums[i]);
                        long double mylg2_bkm = my_log2_bkm(nums[i], bkmd);
                        printf("\n\tref: log2l(%.3LE) = %.20LE\n",
                                nums[i], tgt);
                        printf("\tmy_Lg2(%.3LE)     = %.20LE (- ref: %1.2LE)\n",
                                nums[i], mylg2, (mylg2 - tgt));
                        printf("\tmy_Lg2_bkm(%.3LE) = %.20LE (- ref: %1.2LE)\n",
                                nums[i],
                                mylg2_bkm, mylg2_bkm - tgt);
                }
        }
        */


        return 0;
}