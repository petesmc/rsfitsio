use std::sync::Mutex;

use crate::c_types::{c_int, c_uint};

/*
   Simple Random Number Generators
       - getuniform - uniform deviate [0,1]
       - getnorm    - gaussian (normal) deviate (mean=0, stddev=1)
       - getpoisson - poisson deviate for given expected mean lambda

   This code is adapted from SimpleRNG by John D Cook, which is
   provided in the public domain.

   The original C++ code is found here:
   http://www.johndcook.com/cpp_random_number_generation.html

   This code has been modified in the following ways compared to the
   original.
     1. convert to C from C++
     2. keep only uniform, gaussian and poisson deviates
     3. state variables are module static instead of class variables
     4. provide an srand() equivalent to initialize the state
*/

/*
  These values are not magical, just the default values Marsaglia used.
  Any unit should work.
*/
static M_U: Mutex<c_uint> = Mutex::new(521288629);
static M_V: Mutex<c_uint> = Mutex::new(362436069);

/// Set u and v state variables
fn simplerng_setstate(u: c_uint, v: c_uint) {
    *M_U.lock().unwrap() = u;
    *M_V.lock().unwrap() = v;
}

/// Retrieve u and v state variables
fn simplerng_getstate(u: &mut c_uint, v: &mut c_uint) {
    *u = *M_U.lock().unwrap();
    *v = *M_V.lock().unwrap();
}

/// srand() equivalent to seed the two state variables
pub(crate) fn simplerng_srand(seed: c_uint) {
    fastrand::seed(seed as u64);
}

/// Private routine to get uniform deviate
fn simplerng_getuniform_pr(u: &mut c_uint, v: &mut c_uint) -> f64 {
    /* 0 <= u <= 2^32 */
    let z = simplerng_getuint_pr(u, v);
    /* The magic number is 1/(2^32) and so result is positive and less than 1. */
    (z as f64) * 0.0000000002328306435996595
}

/// Private routine to get unsigned integer
/// Marsaglia multiply-with-carry algorithm (MWC)
fn simplerng_getuint_pr(u: &mut c_uint, v: &mut c_uint) -> c_uint {
    *v = 36969 * ((*v) & 65535) + ((*v) >> 16);
    *u = 18000 * ((*u) & 65535) + ((*u) >> 16);
    ((*v) << 16) + (*u)
}

/// Get uniform deviate [0,1]
pub(crate) fn simplerng_getuniform() -> f64 {
    let r = fastrand::u64(..);

    (r as f64) * (1.0 / ((u64::MAX as f64) + 1.0))
}

/// Get unsigned integer [0, UINT_MAX]
fn simplerng_getuint() -> c_uint {
    /* WARNING: no option for calling rand() here.  Will need to provide
    a scalar to make the uint in the [0,UINT_MAX] range */
    let mut m_u_local = M_U.lock().unwrap();
    let mut m_v_local = M_V.lock().unwrap();
    simplerng_getuint_pr(&mut m_u_local, &mut m_v_local)
}

// Allowed to keep calculation constant between C and Rust version
#[allow(clippy::approx_constant)]
const PI: f64 = 3.141592653589793;

/// Get normal (Gaussian) random sample with mean=0, stddev=1
/// Since you get two deviates for "free" with each calculation, save one of them for later
/// Use Box-Muller algorithm
pub(crate) fn simplerng_getnorm() -> f64 {
    let u1: f64; /* save second value for next call */
    let u2: f64; /* We already saved a value from the last call so use it */
    let r: f64;
    let theta: f64;

    static SAVED: Mutex<bool> = Mutex::new(false);
    static Y: Mutex<f64> = Mutex::new(0.0);

    let mut saved = SAVED.lock().unwrap();
    let mut y = Y.lock().unwrap();

    if !*saved {
        u1 = simplerng_getuniform();
        u2 = simplerng_getuniform();
        r = (-2.0 * u1.ln()).sqrt();
        theta = 2.0 * PI * u2;
        *y = r * theta.cos();
        *saved = true;
        r * theta.sin()
    } else {
        *saved = false;
        *y
    }
}

/// Poisson deviate for expected mean value lambda.
/// lambda should be in the range [0, infinity]
///
/// For small lambda, a simple rejection method is used
/// For large lambda, an approximation is used
pub(crate) fn simplerng_getpoisson(mut lambda: f64) -> c_int {
    if lambda < 0.0 {
        lambda = 0.0;
    }

    if lambda < 15.0 {
        simplerng_poisson_small(lambda)
    } else {
        simplerng_poisson_large(lambda)
    }
}

fn simplerng_poisson_small(lambda: f64) -> c_int {
    /* Algorithm due to Donald Knuth, 1969. */
    let mut p = 1.0;
    let L = (-lambda).exp();
    let mut k = 0;
    loop {
        k += 1;
        p *= simplerng_getuniform();
        if p <= L {
            break;
        };
    }
    k - 1
}

fn simplerng_poisson_large(lambda: f64) -> c_int {
    /* "Rejection method PA" from "The Computer Generation of Poisson Random Variables" by A. C. Atkinson
    Journal of the Royal Statistical Society Series C (Applied Statistics) Vol. 28, No. 1. (1979)
    The article is on pages 29-35. The algorithm given here is on page 32. */

    struct Params {
        beta: f64,
        alpha: f64,
        k: f64,
        old_lambda: f64,
    }

    static PARAMS: Mutex<Params> = Mutex::new(Params {
        beta: 0.0,
        alpha: 0.0,
        k: 0.0,
        old_lambda: -999999.0,
    });

    let mut params = PARAMS.lock().unwrap();

    if lambda != params.old_lambda {
        let c = 0.767 - 3.36 / lambda;
        params.beta = PI / (3.0 * lambda).sqrt();
        params.alpha = params.beta * lambda;
        params.k = (c).ln() - lambda - (params.beta).ln();
        params.old_lambda = lambda;
    }

    loop {
        /* forever */

        let u = simplerng_getuniform();
        let x = (params.alpha - ((1.0 - u) / u).ln()) / params.beta;
        let n = (x + 0.5).floor() as c_int;
        if n < 0 {
            continue;
        }
        let v: f64 = simplerng_getuniform();
        let y: f64 = params.alpha - params.beta * x;
        let temp: f64 = 1.0 + (y).exp();
        let lhs: f64 = y + (v / (temp * temp)).ln();
        let rhs: f64 = params.k + (n as f64) * (lambda).ln() - simplerng_logfactorial(n);
        if lhs <= rhs {
            return n;
        };
    }
}

/// Lookup table for log-gamma function
static LF: [f64; 255] = [
    0.0,
    0.0,
    #[allow(clippy::approx_constant)]
    0.693147180559945,
    1.791759469228055,
    3.178053830347946,
    4.787491742782046,
    6.579251212010101,
    8.525161361065415,
    10.60460290274525,
    12.801827480081469,
    15.104412573075516,
    17.502307845873887,
    19.987214495661885,
    22.55216385312342,
    25.191221182738683,
    27.899271383840894,
    30.671860106080675,
    33.50507345013689,
    36.39544520803305,
    39.339884187199495,
    42.335616460753485,
    45.38013889847691,
    48.47118135183523,
    51.60667556776438,
    54.78472939811232,
    58.00360522298052,
    61.261701761002,
    64.55753862700632,
    67.88974313718153,
    71.257038967168,
    74.65823634883016,
    78.0922235533153,
    81.55795945611503,
    85.05446701758152,
    88.58082754219768,
    92.13617560368708,
    95.7196945421432,
    99.33061245478743,
    102.96819861451381,
    106.63176026064345,
    110.32063971475739,
    114.03421178146169,
    117.77188139974506,
    121.53308151543864,
    125.31727114935688,
    129.12393363912724,
    132.9525750356163,
    136.80272263732635,
    140.67392364823425,
    144.5657439463449,
    148.47776695177302,
    152.40959258449735,
    156.3608363030788,
    160.33112821663093,
    164.32011226319517,
    168.32744544842765,
    172.35279713916282,
    176.39584840699737,
    180.45629141754378,
    184.5338288614495,
    188.6281734236716,
    192.7390472878449,
    196.86618167288998,
    201.00931639928157,
    205.1681994826412,
    209.34258675253682,
    213.53224149456327,
    217.73693411395425,
    221.95644181913036,
    226.19054832372757,
    230.43904356577693,
    234.70172344281826,
    238.97838956183435,
    243.26884900298273,
    247.5729140961869,
    251.8904022097232,
    256.2211355500095,
    260.5649409718632,
    264.9216497985528,
    269.2910976510198,
    273.6731242856937,
    278.0675734403661,
    282.4742926876304,
    286.893133295427,
    291.3239500942703,
    295.7666013507606,
    300.2209486470141,
    304.6868567656687,
    309.1641935801469,
    313.652829949879,
    318.1526396202093,
    322.6634991267262,
    327.1852877037752,
    331.7178871969285,
    336.26118197919845,
    340.81505887079896,
    345.37940706226686,
    349.95411804077025,
    354.5390855194408,
    359.13420536957534,
    363.73937555556347,
    368.3544960724047,
    372.979468885689,
    377.61419787391867,
    382.25858877306,
    386.91254912321756,
    391.5759882173296,
    396.2488170517915,
    400.93094827891576,
    405.6222961611449,
    410.3227765269373,
    415.0323067282496,
    419.7508055995448,
    424.4781934182571,
    429.21439186665157,
    433.95932399501487,
    438.71291418612117,
    443.47508812091894,
    448.2457727453846,
    453.0248962384961,
    457.8123879812781,
    462.6081785268749,
    467.4121995716081,
    472.2243839269805,
    477.0446654925856,
    481.8729792298879,
    486.70926113683936,
    491.553448223298,
    496.4054784872176,
    501.26529089157924,
    506.13282534203483,
    511.00802266523607,
    515.8908245878225,
    520.7811737160442,
    525.679013515995,
    530.5842882944336,
    535.4969431801695,
    540.4169241059977,
    545.344177791155,
    550.2786517242856,
    555.220294146895,
    560.1690540372731,
    565.1248810948744,
    570.0877257251342,
    575.0575390247102,
    580.0342727671308,
    585.0178793888392,
    590.0083119756179,
    595.005524249382,
    600.0094705553274,
    605.0201058494238,
    610.0373856862387,
    615.0612662070849,
    620.0917041284774,
    625.1286567308911,
    630.1720818478102,
    635.2219378550598,
    640.2781836604081,
    645.340778693435,
    650.4096828956552,
    655.4848567108891,
    660.5662610758735,
    665.653857411106,
    670.7476076119127,
    675.8474740397369,
    680.9534195136375,
    686.065407301994,
    691.1834011144108,
    696.307365093814,
    701.4372638087372,
    706.5730622457875,
    711.71472580229,
    716.8622202791034,
    722.0155118736013,
    727.1745671728158,
    732.3393531467393,
    737.5098371417774,
    742.6859868743512,
    747.8677704246434,
    753.0551562304842,
    758.2481130813743,
    763.4466101126402,
    768.650616799717,
    773.8601029525585,
    779.0750387101674,
    784.2953945352457,
    789.521141208959,
    794.7522498258135,
    799.9886917886435,
    805.2304388037031,
    810.4774628758636,
    815.7297363039102,
    820.9872316759379,
    826.2499218648428,
    831.5177800239063,
    836.7907795824699,
    842.0688942417005,
    847.3520979704384,
    852.6403650011331,
    857.9336698258575,
    863.2319871924054,
    868.5352921004646,
    873.8435597978657,
    879.1567657769076,
    884.4748857707518,
    889.7978957498902,
    895.1257719186799,
    900.4584907119453,
    905.7960287916463,
    911.1383630436112,
    916.4854705743288,
    921.8373287078049,
    927.1939149824767,
    932.5552071481862,
    937.9211831632081,
    943.2918211913357,
    948.6670995990198,
    954.0469969525604,
    959.4314920153495,
    964.8205637451659,
    970.2141912915183,
    975.6123539930362,
    981.0150313749084,
    986.4222031463686,
    991.8338491982234,
    997.2499496004278,
    1002.6704845997003,
    1008.0954346171817,
    1013.5247802461362,
    1018.9585022496902,
    1024.3965815586134,
    1029.8389992691355,
    1035.2857366408016,
    1040.7367750943674,
    1046.192096209725,
    1051.6516817238692,
    1057.115513528895,
    1062.58357367003,
    1068.0558443437014,
    1073.5323078956328,
    1079.012946818975,
    1084.4977437524656,
    1089.9866814786224,
    1095.4797429219627,
    1100.976911147256,
    1106.4781693578009,
    1111.983500893733,
    1117.492889230361,
    1123.0063179765261,
    1128.5237708729908,
    1134.045231790853,
    1139.5706847299848,
    1145.100113817496,
    1150.6335033062237,
    1156.1708375732424,
];

fn simplerng_logfactorial(n: c_int) -> f64 {
    if n < 0 {
        return 0.0;
    }
    if n > 254 {
        let x = (n + 1) as f64;
        return (x - 0.5) * (x).ln() - x + 0.5 * (2.0 * PI).ln() + 1.0 / (12.0 * x);
    }
    LF[n as usize]
}
