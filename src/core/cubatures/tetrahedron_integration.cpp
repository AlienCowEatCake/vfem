#include "tetrahedron_integration.h"

#include <cmath>

namespace core { namespace cubatures {

// *************************************************************************************************

namespace {

// Численное интегрирование на тетраэдрах
// https://people.fh-landshut.de/~maurer/numeth/node148.html

namespace tet_integration_2
{
    const std::size_t gauss_num = 4;
    const double gauss_weights[gauss_num] =
    {
        1.0 / 24.0,
        1.0 / 24.0,
        1.0 / 24.0,
        1.0 / 24.0
    };
    const double gauss_a = (5.0 - std::sqrt(5.0)) / 20.0;
    const double gauss_b = (5.0 + 3.0 * std::sqrt(5.0)) / 20.0;
    const double gauss_points_master[gauss_num][4] =
    {
        { 1.0 - gauss_b - 2.0 * gauss_a, gauss_b, gauss_a, gauss_a },
        { 1.0 - gauss_b - 2.0 * gauss_a, gauss_a, gauss_b, gauss_a },
        { 1.0 - gauss_b - 2.0 * gauss_a, gauss_a, gauss_a, gauss_b },
        { 1.0 - 3.0 * gauss_a,           gauss_a, gauss_a, gauss_a }
    };
}

namespace tet_integration_3
{
    const std::size_t gauss_num = 5;
    const double gauss_weights[gauss_num] =
    {
        - 2.0 / 15.0,
        3.0 / 40.0,
        3.0 / 40.0,
        3.0 / 40.0,
        3.0 / 40.0
    };
    const double gauss_a = 1.0 / 4.0;
    const double gauss_b = 1.0 / 2.0;
    const double gauss_c = 1.0 / 6.0;
    const double gauss_points_master[gauss_num][4] =
    {
        { 1.0 - 3.0 * gauss_a,           gauss_a, gauss_a, gauss_a },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_b, gauss_c, gauss_c },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_b, gauss_c },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_c, gauss_b },
        { 1.0 - 3.0 * gauss_c,           gauss_c, gauss_c, gauss_c }
    };
}

namespace tet_integration_4
{
    const std::size_t gauss_num = 11;
    const double gauss_weights[gauss_num] =
    {
        - 74.0 / 5625.0,
        343.0 / 45000.0,
        343.0 / 45000.0,
        343.0 / 45000.0,
        343.0 / 45000.0,
        56.0 / 2250.0,
        56.0 / 2250.0,
        56.0 / 2250.0,
        56.0 / 2250.0,
        56.0 / 2250.0,
        56.0 / 2250.0
    };
    const double gauss_a = 1.0 / 4.0;
    const double gauss_b = 11.0 / 14.0;
    const double gauss_c = 5.0 / 70.0;
    const double gauss_d = (1.0 + std::sqrt(5.0 / 14.0)) / 4.0;
    const double gauss_e = (1.0 - std::sqrt(5.0 / 14.0)) / 4.0;
    const double gauss_points_master[gauss_num][4] =
    {
        { 1.0 - 3.0 * gauss_a,           gauss_a, gauss_a, gauss_a },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_b, gauss_c, gauss_c },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_b, gauss_c },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_c, gauss_b },
        { 1.0 - 3.0 * gauss_c,           gauss_c, gauss_c, gauss_c },
        { 1.0 - gauss_d - 2.0 * gauss_e, gauss_d, gauss_e, gauss_e },
        { 1.0 - gauss_d - 2.0 * gauss_e, gauss_e, gauss_d, gauss_e },
        { 1.0 - gauss_d - 2.0 * gauss_e, gauss_e, gauss_e, gauss_d },
        { 1.0 - gauss_e - 2.0 * gauss_d, gauss_e, gauss_d, gauss_d },
        { 1.0 - gauss_e - 2.0 * gauss_d, gauss_d, gauss_e, gauss_d },
        { 1.0 - gauss_e - 2.0 * gauss_d, gauss_d, gauss_d, gauss_e }
    };
}

} // namespace

// *************************************************************************************************

namespace {

// Интегрирование особо высоких порядков на тетраэдрах
// http://lsec.cc.ac.cn/~tcui/myinfo/paper/quad.pdf
// http://lsec.cc.ac.cn/phg/download.htm

#define tet_points_S4(a)             {(a), (a), (a), (a)}
#define tet_points_S31(a)            {(a), (a), (a), 1.0-3.0*(a)},      {(a), (a), 1.0-3.0*(a), (a)}, \
                                     {(a), 1.0-3.0*(a), (a), (a)},      {1.0-3.0*(a), (a), (a), (a)}
#define tet_points_S22(a)            {(a), (a), 0.5-(a), 0.5-(a)},      {(a), 0.5-(a), (a), 0.5-(a)}, \
                                     {(a), 0.5-(a), 0.5-(a), (a)},      {0.5-(a), (a), 0.5-(a), (a)}, \
                                     {0.5-(a), (a), (a), 0.5-(a)},      {0.5-(a), 0.5-(a), (a), (a)}
#define tet_points_S211(a, b)        {(a), (a), (b), 1.0-(a)-(a)-(b)},  {(a), (a), 1.0-(a)-(a)-(b), (b)}, \
                                     {(a), (b), (a), 1.0-(a)-(a)-(b)},  {(a), (b), 1.0-(a)-(a)-(b), (a)}, \
                                     {(a), 1.0-(a)-(a)-(b), (a), (b)},  {(a), 1.0-(a)-(a)-(b), (b), (a)}, \
                                     {(b), (a), (a), 1.0-(a)-(a)-(b)},  {(b), (a), 1.0-(a)-(a)-(b), (a)}, \
                                     {(b), 1.0-(a)-(a)-(b), (a), (a)},  {1.0-(a)-(a)-(b), (a), (a), (b)}, \
                                     {1.0-(a)-(a)-(b), (a), (b), (a)},  {1.0-(a)-(a)-(b), (b), (a), (a)}
#define tet_points_S0111(p, a, b, c) {(p), (a), (b), (c)},              {(p), (a), (c), (b)}, \
                                     {(p), (b), (a), (c)},              {(p), (b), (c), (a)}, \
                                     {(p), (c), (a), (b)},              {(p), (c), (b), (a)}
#define tet_points_S1111(a, b, c)    tet_points_S0111((a), (b), (c), 1.0-(a)-(b)-(c)), \
                                     tet_points_S0111((b), (a), (c), 1.0-(a)-(b)-(c)), \
                                     tet_points_S0111((c), (a), (b), 1.0-(a)-(b)-(c)), \
                                     tet_points_S0111(1.0-(a)-(b)-(c), (a), (b), (c))
#define tet_weights_S4(w)    (w) / 6.0
#define tet_weights_S31(w)   (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0
#define tet_weights_S22(w)   (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0
#define tet_weights_S211(w)  (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0
#define tet_weights_S0111(w) (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0, (w) / 6.0
#define tet_weights_S1111(w) tet_weights_S0111(w), tet_weights_S0111(w), tet_weights_S0111(w), tet_weights_S0111(w)

namespace tet_integration_5
{
    const std::size_t gauss_num = 14;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S31(0.11268792571801585079918565233328633),
        tet_weights_S31(0.07349304311636194954371020548632750),
        tet_weights_S22(0.04254602077708146643806942812025744)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S31(0.31088591926330060979734573376345783),
        tet_points_S31(0.09273525031089122640232391373703061),
        tet_points_S22(0.04550370412564964949188052627933943)
    };
}

//namespace tet_integration_5_2 // Stroud T3:5-1 p315
//{
//    const std::size_t gauss_num = 15;
//    const double gauss_weights[gauss_num] =
//    {
//        tet_weights_S4(16.0 / 135.0),
//        // (2665 + 14 * std::sqrt(15)) / 37800
//        tet_weights_S31(0.07193708377901862),
//        // (2665 - 14 * std::sqrt(15)) / 37800
//        tet_weights_S31(0.06906820722627239),
//        tet_weights_S22(20.0 / 378.0)
//    };
//    const double gauss_points_master[gauss_num][4] =
//    {
//        tet_points_S4(0.25),
//        // (7 - std::sqrt(15)) / 34, (13 + 3 * std::sqrt(15)) / 34
//        tet_points_S31(0.09197107805272303),
//        // (7 + std::sqrt(15)) / 34, (13 - 3 * std::sqrt(15)) / 34
//        tet_points_S31(0.31979362782962991),
//        // (10 - 2 * std::sqrt(15)) / 40, (10 + 2 * std::sqrt(15)) / 40
//        tet_points_S22(0.05635083268962916)
//    };
//}

namespace tet_integration_6
{
    const std::size_t gauss_num = 24;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S31(0.03992275025816749209969062755747998),
        tet_weights_S31(0.01007721105532064294801323744593686),
        tet_weights_S31(0.05535718154365472209515327785372602),
        tet_weights_S211(27.0 / 560.0)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S31(0.21460287125915202928883921938628499),
        tet_points_S31(0.04067395853461135311557944895641006),
        tet_points_S31(0.32233789014227551034399447076249213),
        tet_points_S211(0.06366100187501752529923552760572698, 0.60300566479164914136743113906093969)
    };
}

namespace tet_integration_7
{
    const std::size_t gauss_num = 35;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S4(0.09548528946413084886057843611722638),
        tet_weights_S31(0.04232958120996702907628617079854674),
        tet_weights_S22(0.03189692783285757993427482408294246),
        tet_weights_S211(0.03720713072833462136961556119148112),
        tet_weights_S211(0.00811077082990334156610343349109654)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S4(0.25),
        tet_points_S31(0.31570114977820279942342999959331149),
        tet_points_S22(0.05048982259839636876305382298656247),
        tet_points_S211(0.18883383102600104773643110385458576, 0.57517163758700002348324157702230752),
        tet_points_S211(0.02126547254148324598883610149981994, 0.81083024109854856111810537984823239)
    };
}

namespace tet_integration_8
{
    const std::size_t gauss_num = 46;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S31(0.00639714777990232132145142033517302),
        tet_weights_S31(0.04019044802096617248816115847981783),
        tet_weights_S31(0.02430797550477032117486910877192260),
        tet_weights_S31(0.05485889241369744046692412399039144),
        tet_weights_S22(0.03571961223409918246495096899661762),
        tet_weights_S211(0.00718319069785253940945110521980376),
        tet_weights_S211(0.01637218194531911754093813975611913)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S31(0.03967542307038990126507132953938949),
        tet_points_S31(0.31448780069809631378416056269714830),
        tet_points_S31(0.10198669306270330000000000000000000),
        tet_points_S31(0.18420369694919151227594641734890918),
        tet_points_S22(0.06343628775453989240514123870189827),
        tet_points_S211(0.02169016206772800480266248262493018, 0.71993192203946593588943495335273478),
        tet_points_S211(0.20448008063679571424133557487274534, 0.58057719012880922417539817139062041)
    };
}

namespace tet_integration_9
{
    const std::size_t gauss_num = 59;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S4(0.05489853459364812686895885032391298),
        tet_weights_S31(0.00421825735654367356185795185819147),
        tet_weights_S31(0.02348412311384798927791501022996111),
        tet_weights_S31(0.00421283454980389148648831814037819),
        tet_weights_S31(0.02994712640542812769203037546126163),
        tet_weights_S22(0.03695441750679136335292416138761121),
        tet_weights_S211(0.00817349224171051348425319650294732),
        tet_weights_S211(0.00987978656102278957913113314297149),
        tet_weights_S211(0.02160718741919244401497646690335203)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S4(0.25000000000000000000000000000000000),
        tet_points_S31(0.03785502061999503609086515586175707),
        tet_points_S31(0.16954439965012220000000000000000000),
        tet_points_S31(0.05484140424416689000000000000000000),
        tet_points_S31(0.32229717190921058836777748445908171),
        tet_points_S22(0.10961777508972033704050355954365052),
        tet_points_S211(0.45915766038590539763886410168178216, 0.08004485927247373376034330857923567),
        tet_points_S211(0.03296694775357210169727386483414899, 0.71879584022434055051132299796383374),
        tet_points_S211(0.18174359672117481549870278661377760, 0.60023700739524674102301240348069459)
    };
}

namespace tet_integration_10
{
    const std::size_t gauss_num = 79;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S4(0.04574189830483037077884770618329337),
        tet_weights_S31(0.01092727610912416907498417206565671),
        tet_weights_S31(0.00055352334192264689534558564012282),
        tet_weights_S31(0.02569337913913269580782688316792080),
        tet_weights_S22(0.00055387649657283109312967562590035),
        tet_weights_S211(0.01044842402938294329072628200105773),
        tet_weights_S211(0.02513844602651287118280517785487423),
        tet_weights_S211(0.01178620679249594711782155323755017),
        tet_weights_S211(0.01332022473886650471019828463616468),
        tet_weights_S211(0.00615987577565961666092767531756180)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S4(0.25000000000000000000000000000000000),
        tet_points_S31(0.11425191803006935688146412277598412),
        tet_points_S31(0.01063790234539248531264164411274776),
        tet_points_S31(0.31274070833535645859816704980806110),
        tet_points_S22(0.01631296303281644000000000000000000),
        tet_points_S211(0.03430622963180452385835196582344460, 0.59830121060139461905983787517050400),
        tet_points_S211(0.12346418534551115945916818783743644, 0.47120066204746310257913700590727081),
        tet_points_S211(0.40991962933181117418479812480531207, 0.16546413290740130923509687990363569),
        tet_points_S211(0.17397243903011716743177479785668929, 0.62916375300275643773181882027844514),
        tet_points_S211(0.03002157005631784150255786784038011, 0.81213056814351208262160080755918730)
    };
}

namespace tet_integration_11
{
    const std::size_t gauss_num = 96;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S31(0.01612698613577620369120244222737879),
        tet_weights_S31(0.00178872341812357138976990346996962),
        tet_weights_S31(0.00847529348343123401863799968389086),
        tet_weights_S31(0.01238021263944669050859562763135516),
        tet_weights_S31(0.02205586697199415746140963638568037),
        tet_weights_S31(0.02295765467664274421265594265203307),
        tet_weights_S22(0.00120553827014535727045055662252294),
        tet_weights_S22(0.02479381575164443454447803302296997),
        tet_weights_S211(0.01203878836480353606935457416590660),
        tet_weights_S211(0.00189370204498242146248858917618493),
        tet_weights_S211(0.01838752922255814184581020943433469),
        tet_weights_S211(0.00375249249801662461193260176157591),
        tet_weights_S211(0.00633289841693951300885921328914879)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S31(0.12460560449278830000000000000000000),
        tet_points_S31(0.02609630765687464746851542316261877),
        tet_points_S31(0.07193883255798884087330011042809557),
        tet_points_S31(0.32611122454203676937273102302894204),
        tet_points_S31(0.29405882789858127213310307732130217),
        tet_points_S31(0.19271399104965490000000000000000000),
        tet_points_S22(0.00047127204692773946587837159205225),
        tet_points_S22(0.10321360207480949336085123341390539),
        tet_points_S211(0.04349989920159741251267172033621503, 0.63045319723555591476353398203997141),
        tet_points_S211(0.01414839289422299290755441603794058, 0.82491678632147090000000000000000000),
        tet_points_S211(0.21646077368258425486341884576246642, 0.52711130286496480000000000000000000),
        tet_points_S211(0.13301884366834711587538262083530116, 0.73318551371398651551736762818473584),
        tet_points_S211(0.44054756810613723082959230959880706, 0.11506799584377921703650823955291194)
    };
}

namespace tet_integration_12
{
    const std::size_t gauss_num = 127;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S4(0.02340581914868067999082580773836836),
        tet_weights_S31(0.00484469946470415656870798306091558),
        tet_weights_S31(0.00079865303812732982185563521014343),
        tet_weights_S31(0.01311872008808756207964488505025527),
        tet_weights_S31(0.02352182961292765917274505054313770),
        tet_weights_S31(0.00210860882494149803857437048649497),
        tet_weights_S31(0.00047839298963616600187228601742259),
        tet_weights_S22(0.00204546234216855322941711800170502),
        tet_weights_S211(0.00334576331671817115245418532677178),
        tet_weights_S211(0.01181044822479275264785338274950585),
        tet_weights_S211(0.00290156990282342152841364375092118),
        tet_weights_S211(0.00949250645501753676094846901252898),
        tet_weights_S211(0.02094018358085748583183796760479700),
        tet_weights_S211(0.00171435866337409051521874943702732),
        tet_weights_S1111(0.00759915954173370886076474450830409)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S4(0.25000000000000000000000000000000000),
        tet_points_S31(0.19318721110347230000000000000000000),
        tet_points_S31(0.01811701371436566878506928822499717),
        tet_points_S31(0.10700751831426066518406159227423033),
        tet_points_S31(0.29936173715970702940603127680004538),
        tet_points_S31(0.33333033333333333042835213613025030),
        tet_points_S31(0.16575369007421640000000000000000000),
        tet_points_S22(0.04009986052352575650366980228640728),
        tet_points_S211(0.01951844463761131301132122485607343, 0.59982639757597731668263005976738196),
        tet_points_S211(0.24970741896308715787490891769354198, 0.47400425629911050000000000000000000),
        tet_points_S211(0.07674205857869954726322831328843659, 0.83056291375422969598432041821082569),
        tet_points_S211(0.43011409627915217536723647418133112, 0.02265922072588833582931396831630072),
        tet_points_S211(0.12197854304894211937147375564906792, 0.47765370899783134571567376444973682),
        tet_points_S211(0.01480482319031682427540691439704854, 0.81083799468092699988474915243749073),
        tet_points_S1111(0.65250697573013212016385330106711095, 0.22646235632397177636617160407210034, 0.02251830769546778956654013747639605)
    };
}

namespace tet_integration_13
{
    const std::size_t gauss_num = 149;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S4(0.02191579945212728678229670892998658),
        tet_weights_S31(0.00809592740005652573580359966615063),
        tet_weights_S31(0.00130319185047278813746994806952476),
        tet_weights_S31(0.01996610676014222116016391561580003),
        tet_weights_S31(0.02125705756007566772097136088386650),
        tet_weights_S22(0.00077331890737182713690269661719116),
        tet_weights_S22(0.01755491389570430512641028370006205),
        tet_weights_S211(0.00213830361001659899343287397434178),
        tet_weights_S211(0.00256560169283338620814651902766716),
        tet_weights_S211(0.00338953948455728203040932651810398),
        tet_weights_S211(0.01135828330503278417235563981454793),
        tet_weights_S211(0.01103203882197761043040360052454856),
        tet_weights_S211(0.00457602573785952356043458354199517),
        tet_weights_S211(0.00827343104220868129752243222682095),
        tet_weights_S211(0.00586641165391940007076312979369247),
        tet_weights_S1111(0.00313458521939849614410720196518793)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S4(0.25000000000000000000000000000000000),
        tet_points_S31(0.09935339765028269917868020572165369),
        tet_points_S31(0.02361873260499568532036302265004401),
        tet_points_S31(0.30089166537572662790706731844610997),
        tet_points_S31(0.18156624280757148139366685840064601),
        tet_points_S22(0.00428160639152879988718710754508354),
        tet_points_S22(0.12290357421888442998582785890620434),
        tet_points_S211(0.28318219770202728236417353077594322, 0.43037955664247500440987356786807501),
        tet_points_S211(0.02239485904524970717572425710098278, 0.83488749018470024820940398932904512),
        tet_points_S211(0.02191788402113435132324662419880111, 0.67691762094326571059673391273529166),
        tet_points_S211(0.21481417044274656673534260788169227, 0.52280311286258745560867693994038579),
        tet_points_S211(0.08000490008644308882018405418010744, 0.24689045570275147370034631014113188),
        tet_points_S211(0.11579466150271899371721034492503850, 0.74997281767443310000000000000000000),
        tet_points_S211(0.39129315347000474438672195978809687, 0.18835457382799180000000000000000000),
        tet_points_S211(0.45315745821242834581317282468854978, 0.02202033169457796534173826092007299),
        tet_points_S1111(0.27324999892429634023602493512400674, 0.60775441245653315696274741541102470, 0.00561877924700169073874366184065955)
    };
}

namespace tet_integration_14
{
    const std::size_t gauss_num = 236;
    const double gauss_weights[gauss_num] =
    {
        tet_weights_S31(0.00406511366527076704362088368356360),
        tet_weights_S31(0.00221453853344557814375995695000715),
        tet_weights_S31(0.00581343826788845054953733388214554),
        tet_weights_S31(0.01962554338583572159756233339617148),
        tet_weights_S31(0.00038757379059082143645387212483937),
        tet_weights_S211(0.01164297197217703698552134010055516),
        tet_weights_S211(0.00528904298828171313177368830528561),
        tet_weights_S211(0.00183108541636005593766978234880692),
        tet_weights_S211(0.00824964737721464520674496691736603),
        tet_weights_S1111(0.00300992453470824513768887482089866),
        tet_weights_S1111(0.00080471656173675346362618087603116),
        tet_weights_S1111(0.00298504125884930711876556928839215),
        tet_weights_S1111(0.00568960024187607669633614778119730),
        tet_weights_S1111(0.00415908658785457156700139801826135),
        tet_weights_S1111(0.00072823892045727243561364297456536),
        tet_weights_S1111(0.00543265007699582482162423406519264)
    };
    const double gauss_points_master[gauss_num][4] =
    {
        tet_points_S31(0.32725336252384856390930966926852893),
        tet_points_S31(0.04476130446668508088379420964788419),
        tet_points_S31(0.08614033110243635365372087402988575),
        tet_points_S31(0.20876264250043229682653570839761758),
        tet_points_S31(0.01410497380292096006358791521029282),
        tet_points_S211(0.10216532418077681234766925269825839, 0.57394636759433382028140028934601068),
        tet_points_S211(0.40757005166001071572132956513017833, 0.09222787013902013000000000000000000),
        tet_points_S211(0.01566400074028035855575867095780840, 0.70128109595894403271399676732084261),
        tet_points_S211(0.22549635625250290537807241542011034, 0.47690639744208871158605833541070112),
        tet_points_S1111(0.39059842812814580000000000000000000, 0.20135905441239221681230773272350923, 0.01611228807103002985780269315483708),
        tet_points_S1111(0.10613506799890214555561390298480794, 0.03273581868172692849440040779126601, 0.00359790765372716669079715233859245),
        tet_points_S1111(0.56363837316977438968968166306485017, 0.23029207223006574545025268741356515, 0.19071993417435518627124877906378985),
        tet_points_S1111(0.36762550953258608440922067759911669, 0.20788513802300449507171021252507348, 0.33121048851934490000000000000000000),
        tet_points_S1111(0.71923236898172952950234018407969909, 0.17632791180193297621579930336369727, 0.02076023625713100907549734406116442),
        tet_points_S1111(0.52782499521529872984092400758172763, 0.43728908922034181655262387608419181, 0.00922016518566419494631775549492202),
        tet_points_S1111(0.54836745449481907289949105056077457, 0.34478155061716412287036718709203314, 0.08672172833222153946294387400858277)
    };
}

} // namespace

// *************************************************************************************************

void tetrahedron_integration::set(std::size_t num, const double weights[], const double points[][4])
{
    m_gauss_num = num;
    m_gauss_weights.resize(m_gauss_num);
    m_gauss_points_master.resize(m_gauss_num, 4);
    for(std::size_t i = 0; i < m_gauss_num; i++)
    {
        m_gauss_weights[i] = weights[i];
        for(std::size_t j = 0; j < 4; j++)
            m_gauss_points_master[i][j] = points[i][j];
    }
}

tetrahedron_integration::tetrahedron_integration(std::size_t order)
{
    init(order);
}

void tetrahedron_integration::init(std::size_t order)
{
    if(order < 2)
        order = 2;
    if(order > 14)
        order = 14;

    switch(order)
    {
    case 2:
    {
        namespace curr = tet_integration_2;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 3:
    {
        namespace curr = tet_integration_3;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 4:
    {
        namespace curr = tet_integration_4;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 5:
    {
        namespace curr = tet_integration_5;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 6:
    {
        namespace curr = tet_integration_6;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 7:
    {
        namespace curr = tet_integration_7;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 8:
    {
        namespace curr = tet_integration_8;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 9:
    {
        namespace curr = tet_integration_9;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 10:
    {
        namespace curr = tet_integration_10;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 11:
    {
        namespace curr = tet_integration_11;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 12:
    {
        namespace curr = tet_integration_12;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 13:
    {
        namespace curr = tet_integration_13;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 14:
    {
        namespace curr = tet_integration_14;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    default:
    {
        namespace curr = tet_integration_8;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    }
}

// *************************************************************************************************

}} // namespace core::cubatures
