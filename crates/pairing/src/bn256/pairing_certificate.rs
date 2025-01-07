use ff::{Field, PrimeField};
use num_traits::identities::Zero;
use rand::Rand;
use crate::{bn256::{Bn256, G1Affine, G2Affine}, compact_bn256, CurveAffine, Engine, CurveProjective};
use crate::bn256::*;
use crate::bn256::fq::*;
// use crate::bn256::FROBENIUS_COEFF_FQ6_C1;
// use crate::bn256::XI_TO_Q_MINUS_1_OVER_2;

use super::{Fq, Fq2, Fq12, Fq6, FqRepr};
use num_bigint::BigUint;

// power = (|Fq12| - 1)/3
const FQ12_MULT_ORDER_DIV_3 : [u64; 48] = [
    16953508628923568224, 13877003111090065395, 13157049836381475239, 11415053979013103756, 8384535856576305654, 16858761491289863955, 5606204346564449596, 
    11477821246166416167, 12793023930672435534, 6172907199131808173, 6931315615918879492, 2145656818538089052, 4588009837644227054, 14863053780587399191, 
    12417100097682464471, 6028805810388919583, 1386536852145720983, 9857992206758174354, 10109427491521724439, 7101431254775646169, 12194607637311870232, 
    4123599041920339461, 4361107823611899910, 17615083610772831412, 10982264379420315624, 6419188508526177463, 7498391412603358138, 12991551331951386245, 
    11052497498275195037, 11785136137640826697, 636005296143811820, 13105056392368402529, 3966494654182370702, 17822188952467346335, 17127093937586648918, 
    6576355335316381401, 16739546170660304941, 11212174993167722634, 1865707425941064056, 1346562786881031727, 12688767027439644562, 1866376690783873247, 
    16955186452462703345, 14122185803590355033, 15354689103052496694, 7224017006226930792, 8327800425992802905, 12799192649
];

const FQ12_MULT_ORDER : [u64; 48] = [
    13967037739351601440, 4737521185851092955, 2577661361725322487, 15798417863329759654, 6706863496019365347, 13682796326450488634, 16818613039693348790, 
    15986719664789696885, 1485583644598203371, 71977523685872905, 2347202774047086861, 6436970455614267157, 13764029512932681162, 7695673194343094341, 
    357812145628290183, 18086417431166758751, 4159610556437162949, 11127232546564971446, 11881538400855621702, 2857549690617386892, 18137078838226059081, 
    12370797125761018384, 13083323470835699730, 15951762684899391004, 14500049064551395258, 810821451868980774, 4048430164100522799, 2081165848435055504, 
    14710748421116033497, 16908664339212928476, 1908015888431435461, 2421681029686104355, 11899483962547112108, 16573078709982935773, 14487793665340843524, 
    1282321932239592589, 13325150364561811592, 15189780905793616288, 5597122277823192169, 4039688360643095181, 1172812934899830454, 5599130072351619743, 
    13972071209969006803, 5473069263351961869, 9170579161738386852, 3225306944971240762, 6536657204268857100, 38397577948
];

const FQ6_MULT_ORDER_DIV_27 : [u64; 24] = [
    8168956762542115120, 12803991695383604623, 16702320209049413708, 17652318445546162000, 16303848130085388482, 7056390770558977116, 7260638814434544408, 
    6439226801849030222, 17544709870344963018, 11956284243821547882, 5575926253511864029, 7595604974183140837, 11723832411009852653, 12634391244340473392, 
    7019961605730375294, 12764185526556718837, 16666459039236151248, 14472658015929537488, 8140829297918308879, 14928924455053179571, 905557763412333626, 
    10210676808216562665, 9677457426615697220, 31170804581356
];

const FQ6_MULT_ORDER : [u64; 24] = [
    17647647777832040464, 13666382448585395744, 8240787875304931350, 15443996187007583624, 15928785816985801871, 6055110067996865995, 11569807252637182866, 
    7838426986537851460, 12538564656575211095, 9225025330119415367, 2976056255143915872, 2167149492139734831, 2948825844203644170, 9087170270420852513, 
    5071522617624616796, 12591615890259479521, 7272536290346844930, 3380140882196928264, 16888206232989271978, 15699334738535264492, 6003315538423456307, 
    17433856789913469332, 3036933486690102330, 841611723696626
];

const R_PRIME : [u64; 44] = [
    11608900242075462817, 11957288855914354513, 3116787106681700788, 4056127484362567420, 8118514525076557193, 69240188769302606, 9133271292293368251, 
    2039292154691536571, 14287587744049093151, 3554310640973855491, 12261137717040107045, 2251832533098411731, 7336272344092593252, 14765170906362374105, 
    15300949932261141727, 14106183133602261101, 12737401427063335353, 171177352230438000, 13049041856653590179, 18305798765839010998, 16104228988624323963, 
    1739542763929064922, 5601136942292264253, 15574900073644834590, 18274434499964322805, 17916135654498903305, 18322269162143131602, 4007644745573738693, 
    13947846285796705775, 16368008805899974858, 472327754523454952, 5885354980194466143, 2546607790591198545, 2891428736624824818, 17184353915837753597, 
    14744553112369412456, 9373784872161472250, 16200138677860926851, 11959944834529249218, 9547417553780286347, 5489829204756096954, 17399753614172681083, 
    8755379657557813582, 182295210734
];

const M_PRIME : [u64; 48] = [
    10887798433163368087, 6442620111203252320, 5925498864076479409, 16024781233117228933, 12017805302040997602, 16791140354641337449, 12674593397576570157, 
    7215868095803522290, 15636602336049209438, 16730884744177924545, 18268533137212137744, 15527882843810836270, 7081124931869209512, 9516729772367226167, 
    10183752447065015983, 15502457473202806646, 13085989172809235752, 108420017516737598, 1284723564304771983, 8061989650089060148, 7954858051136240117, 
    17888755599543366389, 11649540998901049035, 5945514820183971169, 11922621809996018633, 1933202901127852564, 7078564943103599682, 7718100670145644596, 
    1364989304968683823, 15395144788315683924, 4577511551603766591, 3201951561493279966, 4700210007183877234, 14860867369575060718, 18421873388246629180, 
    6029614175412069092, 2152860897591860129, 10387654856425849196, 15895727615005279996, 14248830315588959169, 1132513893817127473, 1035096278108215019, 
    9772924880967139577, 11715034521440154224, 436422959846410301, 9762853756639923596, 6159059254438582067, 11919040216
];

// write |Fq12| - 1 = 3^t * s, where s is not divisible by t; then TONELLI_SHANKS_EXP = (s + 1) / 3
const TONELLI_SHANKS_EXP: [u64; 48] = [
    8143247905175134667, 1880388565129969208, 14834765754973409599, 5205268981295554261, 15341218721414312637, 9506164238870890183, 10455828350081767549, 
    16822210333896135368, 473815701136016130, 7060753627267678679, 14604183006067387534, 13060510897000724805, 7685266468461084993, 5332972677650157796, 
    7975232774388427120, 6372203795324995708, 17131671803514241162, 1731536309413973243, 16088316340253385615, 13244058098342856550, 12749481517188288863, 
    4252002351265838857, 3577586229339246592, 5434899708397766397, 15437430888927053747, 8436300644186696179, 15991611300293446122, 14828636180735258154, 
    11340755654726963736, 10684677675677188923, 17103874338477133415, 1851797945918055768, 16544013052711541092, 18423612402552432901, 2683975042915381620, 
    6392483407359346146, 16333876291332592300, 9980244149077831305, 16466206118332233438, 5515574643576201653, 9351719999469030206, 5534827010387418006, 
    2677608099021902155, 7355171353358735970, 15599372545357875268, 7099683620123053590, 3724500770168168925, 474044172
];

const S: [u64; 48] = [
    5982999641815852384, 5641165695389907625, 7610809117501125565, 15615806943886662785, 9130168016823834679, 10071748642903118935, 12920740976535751032, 
    13573142854269302873, 1421447103408048392, 2735516808093484421, 6919060870783059371, 2288044543583071185, 4609055331673703365, 15998918032950473389, 
    5478954249455729744, 669867312265435509, 14501527263123620255, 5194608928241919731, 11371460873341053613, 2838686147609466420, 1354956404145763359, 
    12756007053797516573, 10732758688017739776, 16304699125193299191, 9418804519362058009, 6862157858850536923, 11081345753461235135, 7592420394786671232, 
    15575522890471339594, 13607288953322015154, 14418134868012297014, 5555393837754167306, 12738551010715520044, 18377349060238195473, 8051925128746144862, 
    730706148368486822, 12108140726578673669, 11493988373523942301, 12505130207577597083, 16546723930728604961, 9608415924697539002, 16604481031162254019, 
    8032824297065706465, 3618769986366656294, 9904629488654522573, 2852306786659609156, 11173502310504506776, 1422132516
];

const LAMBDA : [u64; 12] = [
    6338883757087263829, 18082362518251613329, 12769575203635315010, 16977236642651234352, 8686951213117475924, 13764436133107921006, 12359884445821460169, 
    15944377745085967733, 13997781615140285916, 7499777415441582065, 3683594101071511471, 124599342199167405
];


#[test]
fn check_cubic_non_residue() {
    // logic is the following: if b = a^3, then b^power = a^(|Fq12| - 1) = 1
    // so we just assemble ranodm element and assumr it's pow is not equal to one
    let elem: Fq12 = FQ12_CUBIC_NON_RESIDUE;
    assert_eq!(elem.pow(&FQ12_MULT_ORDER), Fq12::one());
    
    let pow = elem.pow(&FQ12_MULT_ORDER_DIV_3);
    assert_ne!(pow, Fq12::one()); 
}


#[test]
fn check_root_27_of_unity() {
    // the idea is to start with random element of Fp6, and raise it to pow = (|Fp6| - 1)/27
    // let the new elemetn to be a, it's order is 1, 3, 9 or 27, so we check that a^9 != 1 and assert that a^27 == 1
    let candidate_root_of_unity = ROOT_27_OF_UNITY;
    assert_eq!(candidate_root_of_unity.pow(&[27]), Fq6::one());
    assert_ne!(candidate_root_of_unity.pow(&[9]), Fq6::one());
}


#[derive(Debug)]
pub struct Certificate {
    pub c: Fq12,
    pub root_27_of_unity_power: usize
}


fn cube_elem(elem: &Fq12) -> Fq12 {
    let mut res = elem.clone();
    res.square();
    res.mul_assign(elem);

    res
}


// given element a of order 3^t find t by succesive cubing
fn get_cubic_ord(mut a: Fq12) -> BigUint {
    let mut result = BigUint::ZERO;
    while a != Fq12::one() {
        a = cube_elem(&a);
        result += 1u64;
    }

    result
}


pub fn construct_certificate(miller_loop_f: Fq12) -> Certificate {
    // Input: Output of a Miller loop f and fixed 27-th root of unity w
    // Output: (c, w^i), such that c^λ = f · w^i
    // 1) for i in 0, 1, 2, find the only i, such that (f * w^i)^power == 1
    let mut correct_w_power = 3;
    let mut f = miller_loop_f;

    for i in 0..3 {
        if f.pow(FQ12_MULT_ORDER_DIV_3) == Fq12::one() {
            correct_w_power = i;
            break;
        } else {
            f.c0.mul_assign(&ROOT_27_OF_UNITY);
            f.c1.mul_assign(&ROOT_27_OF_UNITY); 
        }
    }
    assert!(correct_w_power < 3);

    // c = f^r_prime
    let mut c = f.pow(&R_PRIME);
    
    // c = c^m_prime
    c = c.pow(&M_PRIME);
    
    let c_inv = c.inverse().unwrap();
    let cubic_non_residue_pow = FQ12_CUBIC_NON_RESIDUE.pow(&S);

    // extract cubic root of c by using modified Tonelli-Shanks
    let mut x = c.pow(&TONELLI_SHANKS_EXP);
    
    // compute a = x^3 / c
    let mut a = cube_elem(&x);
    a.mul_assign(&c_inv);
    let mut t : BigUint = get_cubic_ord(a);

    while !t.is_zero() {
        x.mul_assign(&cubic_non_residue_pow);
        
        let mut a = cube_elem(&x);
        a.mul_assign(&c_inv);
        t = get_cubic_ord(a);
    }

    let certificate = Certificate {
        c: x,
        root_27_of_unity_power: correct_w_power
    };

    certificate
}

pub fn validate_ceritificate(miller_loop_f: &Fq12, certificate: &Certificate) -> bool {
    // c^λ =? f · w^i
    let lhs = certificate.c.pow(&LAMBDA);
    let w = ROOT_27_OF_UNITY.pow(&[certificate.root_27_of_unity_power as u64]);
    let mut rhs = miller_loop_f.clone();
    rhs.c0.mul_assign(&w);
    rhs.c1.mul_assign(&w);

    lhs == rhs
}


pub fn prepare_all_line_functions(q: G2Affine) -> Vec<(Fq2, Fq2)> {
    fn line_double(t: G2Affine) -> (Fq2, Fq2) {
        let mut alpha = t.x;
        let mut mu = t.y;
    
        // alpha = 3 * x^2 / 2 * y
        alpha.square();
        let mut tmp = alpha.clone();
        tmp.double();
        alpha.add_assign(&tmp);

        let mut tmp = t.y;
        tmp.double();
        let tmp_inv = tmp.inverse().unwrap();
        alpha.mul_assign(&tmp_inv);

        let mut tmp = t.x;
        tmp.mul_assign(&alpha);
        mu.sub_assign(&tmp);

        (alpha, mu)
    }

    fn line_add(t: G2Affine, p: G2Affine) -> (Fq2, Fq2) {
        let x1 = t.x;
        let y1 = t.y;
        let x2 = p.x;
        let y2 = p.y;
    
        let mut alpha = y2;
        let mut mu = y1; 

        // alpha = (y2 - y1) / (x2 - x1)
        alpha.sub_assign(&y1);
        let mut tmp = x2; 
        tmp.sub_assign(&x1);
        let tmp_inv = tmp.inverse().unwrap();
        alpha.mul_assign(&tmp_inv);

        let mut tmp = x1; 
        tmp.mul_assign(&alpha);
        mu.sub_assign(&tmp);
    
        (alpha, mu)
    }

    let mut l = vec![];
    let mut t = q.into_projective();
    let mut q_negated = q.clone();
    q_negated.negate();

    for i in (1..SIX_U_PLUS_2_NAF.len()).rev()  {
        let (alpha, mu) = line_double(t.into_affine());
        t.double();
        l.push((alpha, mu));

        let bit = SIX_U_PLUS_2_NAF[i-1];

        if bit != 0 {
            let q_t = if bit == 1 { q } else { q_negated };
            let (alpha, mu) = line_add(t.into_affine(), q_t);
            t.add_assign_mixed(&q_t);
            l.push((alpha, mu));
        }
    }

    // Frobenius map calculations for BN256
    let mut pi_1_q_x = q.x; 
    let mut pi_1_q_y = q.y; 
    pi_1_q_x.conjugate();
    pi_1_q_x.mul_assign(&FROBENIUS_COEFF_FQ6_C1[1]);
    pi_1_q_y.conjugate();
    pi_1_q_y.mul_assign(&XI_TO_Q_MINUS_1_OVER_2);
    let pi_1_q = G2Affine::from_xy_checked(pi_1_q_x, pi_1_q_y).unwrap();

    let mut pi_2_q_x = q.x; 
    pi_2_q_x.mul_assign(&FROBENIUS_COEFF_FQ6_C1[2]);
    let pi_2_q = G2Affine::from_xy_checked(pi_2_q_x, q.y).unwrap();

    let (alpha, mu) = line_add(t.into_affine(), pi_1_q);
    t.add_assign_mixed(&pi_1_q);
    l.push((alpha.into(), mu.into()));

    let (alpha, mu) = line_add(t.into_affine(), pi_2_q);
    t.add_assign_mixed(&pi_2_q);
    l.push((alpha.into(), mu.into()));

    l
}

pub fn prepare_g1_point(p: G1Affine) -> G1Affine {
    // we "prepare" p by recomputing (x, y) -> (x', y') where: x' = - p.x / p.y; y' = -1 /p.y 
    // it is required to enforce that in c0c3c4, c0 = 1
    let mut y_new = p.y.inverse().unwrap();
    y_new.negate();
    let mut x_new = p.x;
    x_new.mul_assign(&y_new);
    
    G1Affine::from_xy_unchecked(x_new, y_new)
}


pub fn miller_loop_with_prepared_lines(
    eval_points: &[G1Affine],
    lines: &[Vec<(Fq2, Fq2)>],
) -> Fq12 {
    assert_eq!(eval_points.len(), lines.len());

    fn line_evaluation(alpha: &Fq2, mu: &Fq2, p: &G1Affine) -> (Fq2, Fq2, Fq2) {
        // previously: c0 = p.y; c3 = - lambda * p.x; c4 = -mu; 
        // now we have: p = (x', y') where: x' = - p.x / p.y; y' = -1 /p.y 
        // and compute c0 = 1, c3 = - lambda * p.x / p.y = lambda * p.x, c4 = - mu / p.y = mu * p.y
        
        // previously:
        // let mut c3 = *alpha;
        // c3.negate();
        // c3.c0.mul_assign(&p.x);
        // c3.c1.mul_assign(&p.x);
        // let mut c4 = *mu;
        // c4.negate();

        let mut c3 = *alpha;
        c3.c0.mul_assign(&p.x);
        c3.c1.mul_assign(&p.x);
        let mut c4 = *mu;
        c4.c0.mul_assign(&p.y);
        c4.c1.mul_assign(&p.y);
        let c0 = Fq2::one();

        (c0, c3, c4)
    }
    
    let mut f = Fq12::one();

    let mut lc = 0; 
    for i in (1..SIX_U_PLUS_2_NAF.len()).rev() {
        if i != SIX_U_PLUS_2_NAF.len() - 1 {
            f.square();
        }
        let x = SIX_U_PLUS_2_NAF[i - 1];
        for  (P, L) in eval_points.iter().zip(lines.iter()) {
            let (alpha, mu) = L[lc];
            let (c0, c1, c2) = line_evaluation(&alpha, &mu, P);
            f.mul_by_034(&c0, &c1, &c2);

            if x * x == 1 {
                let (alpha, bias) = L[lc + 1];
                let (c0, c1, c2) = line_evaluation(&alpha, &bias, P);
                f.mul_by_034(&c0, &c1, &c2);
            }
        }

        if x == 0 {
            lc += 1;
        } else {
            lc += 2;
        }
    }

    // Frobenius map part: p - p^2 + p^3 steps
    // this part runs through each eval point and applies
    // three additional line evaluations with special conditions.
    // Todo_O_O need to ckeck if in circuits 3 frob line eval or 2 ???
    for (P,  L) in eval_points.iter().zip(lines.iter()) {
        for k in 0..2 {
            let (alpha, mu) = L[lc + k];
            
            let (c0, c1, c2) = line_evaluation(&alpha, &mu, P);
            f.mul_by_034(&c0, &c1, &c2);
        }
    }

    lc += 2;

    // Check we consumed all lines
    if !lines.is_empty() {
        assert_eq!(lc, lines[0].len());
    }

    f
}


#[test]
fn test_certificate_construction() {
    let mut rng = rand::thread_rng();
    let g1 = G1Affine::rand(&mut rng);
    let g2 = G2Affine::rand(&mut rng);

    let g1_prepared = g1.prepare();
    let g2_prepared = g2.prepare();
    
    let mut g2_negated = g2.clone();
    g2_negated.negate();
    let g2_negated_prepared = g2_negated.prepare();

    let miller_loop_result = Bn256::miller_loop(&[(&g1_prepared, &g2_prepared), (&g1_prepared, &g2_negated_prepared)]);
    let final_exp = Bn256::final_exponentiation(&miller_loop_result).unwrap();
    assert_eq!(final_exp, Fq12::one());

    let certificate = construct_certificate(miller_loop_result);
    assert!(validate_ceritificate(&miller_loop_result, &certificate));

    let miller_loop_result = Bn256::miller_loop(&[(&g1_prepared, &g2_prepared)]);
    let final_exp = Bn256::final_exponentiation(&miller_loop_result).unwrap();
    assert_ne!(final_exp, Fq12::one());

    let certificate = construct_certificate(miller_loop_result);
    assert!(!validate_ceritificate(&miller_loop_result, &certificate));
}

#[test]
fn test_bn_equivalence() {
    let mut rng = rand::thread_rng();
    let g1 = G1Affine::rand(&mut rng);
    let g2 = G2Affine::rand(&mut rng);
    let miller_loop_result = Bn256::miller_loop(&[(&g1.prepare(), &g2.prepare())]);
    let lhs = Bn256::final_exponentiation(&miller_loop_result).unwrap();

    let g1 = compact_bn256::G1Affine::from_xy_checked(g1.x, g1.y).unwrap();
    let g2 = compact_bn256::G2Affine::from_xy_checked(g2.x, g2.y).unwrap();
    let miller_loop_result = compact_bn256::Bn256::miller_loop(&[(&g1.prepare(), &g2.prepare())]);
    let rhs = compact_bn256::Bn256::final_exponentiation(&miller_loop_result).unwrap();

    assert_eq!(lhs, rhs)
}


#[test]
fn test_precomputed_lines() {
    let thread_rng = rand::thread_rng();
    let mut rng = thread_rng;
    let p = G1Affine::rand(&mut rng);
    let q = G2Affine::rand(&mut rng);

    let res = Bn256::miller_loop(&[(&p.prepare(), &q.prepare())]);
    let lhs = Bn256::final_exponentiation(&res);

    let lines = prepare_all_line_functions(q);
    let p_prepared = prepare_g1_point(p);
   
    let miller_loop_f = miller_loop_with_prepared_lines(&[p_prepared], &[lines]);
    let rhs = Bn256::final_exponentiation(&miller_loop_f);

    assert_eq!(lhs, rhs);

    let certificate = construct_certificate(miller_loop_f);
    assert!(!validate_ceritificate(&miller_loop_f, &certificate));
}

