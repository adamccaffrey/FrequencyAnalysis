static const float hann[2048] =
{
	0.00000000,
	0.00000236,
	0.00000942,
	0.00002120,
	0.00003769,
	0.00005888,
	0.00008479,
	0.00011541,
	0.00015074,
	0.00019077,
	0.00023552,
	0.00028498,
	0.00033914,
	0.00039801,
	0.00046159,
	0.00052987,
	0.00060286,
	0.00068056,
	0.00076295,
	0.00085006,
	0.00094186,
	0.00103837,
	0.00113958,
	0.00124549,
	0.00135610,
	0.00147140,
	0.00159140,
	0.00171610,
	0.00184549,
	0.00197958,
	0.00211836,
	0.00226183,
	0.00240999,
	0.00256283,
	0.00272037,
	0.00288259,
	0.00304949,
	0.00322107,
	0.00339734,
	0.00357828,
	0.00376390,
	0.00395420,
	0.00414917,
	0.00434881,
	0.00455312,
	0.00476210,
	0.00497574,
	0.00519405,
	0.00541702,
	0.00564465,
	0.00587694,
	0.00611389,
	0.00635548,
	0.00660173,
	0.00685263,
	0.00710817,
	0.00736836,
	0.00763318,
	0.00790265,
	0.00817675,
	0.00845549,
	0.00873886,
	0.00902685,
	0.00931947,
	0.00961672,
	0.00991858,
	0.01022507,
	0.01053616,
	0.01085187,
	0.01117219,
	0.01149711,
	0.01182664,
	0.01216076,
	0.01249948,
	0.01284279,
	0.01319070,
	0.01354319,
	0.01390026,
	0.01426191,
	0.01462814,
	0.01499894,
	0.01537432,
	0.01575425,
	0.01613875,
	0.01652781,
	0.01692143,
	0.01731959,
	0.01772230,
	0.01812956,
	0.01854136,
	0.01895769,
	0.01937855,
	0.01980395,
	0.02023386,
	0.02066830,
	0.02110725,
	0.02155072,
	0.02199869,
	0.02245117,
	0.02290815,
	0.02336962,
	0.02383558,
	0.02430603,
	0.02478096,
	0.02526036,
	0.02574424,
	0.02623259,
	0.02672540,
	0.02722267,
	0.02772440,
	0.02823057,
	0.02874119,
	0.02925625,
	0.02977574,
	0.03029967,
	0.03082802,
	0.03136079,
	0.03189797,
	0.03243957,
	0.03298557,
	0.03353597,
	0.03409077,
	0.03464995,
	0.03521352,
	0.03578147,
	0.03635379,
	0.03693048,
	0.03751154,
	0.03809695,
	0.03868671,
	0.03928082,
	0.03987927,
	0.04048206,
	0.04108917,
	0.04170061,
	0.04231636,
	0.04293643,
	0.04356081,
	0.04418948,
	0.04482245,
	0.04545971,
	0.04610125,
	0.04674706,
	0.04739715,
	0.04805150,
	0.04871011,
	0.04937297,
	0.05004008,
	0.05071142,
	0.05138700,
	0.05206681,
	0.05275083,
	0.05343907,
	0.05413152,
	0.05482817,
	0.05552901,
	0.05623404,
	0.05694325,
	0.05765664,
	0.05837419,
	0.05909590,
	0.05982177,
	0.06055178,
	0.06128594,
	0.06202423,
	0.06276664,
	0.06351317,
	0.06426382,
	0.06501857,
	0.06577742,
	0.06654036,
	0.06730739,
	0.06807849,
	0.06885366,
	0.06963289,
	0.07041618,
	0.07120352,
	0.07199489,
	0.07279030,
	0.07358973,
	0.07439318,
	0.07520064,
	0.07601210,
	0.07682756,
	0.07764700,
	0.07847042,
	0.07929782,
	0.08012918,
	0.08096449,
	0.08180375,
	0.08264695,
	0.08349409,
	0.08434515,
	0.08520012,
	0.08605900,
	0.08692179,
	0.08778846,
	0.08865902,
	0.08953345,
	0.09041175,
	0.09129391,
	0.09217992,
	0.09306977,
	0.09396346,
	0.09486097,
	0.09576230,
	0.09666744,
	0.09757638,
	0.09848911,
	0.09940562,
	0.10032590,
	0.10124996,
	0.10217776,
	0.10310932,
	0.10404462,
	0.10498364,
	0.10592639,
	0.10687285,
	0.10782302,
	0.10877688,
	0.10973442,
	0.11069565,
	0.11166054,
	0.11262909,
	0.11360129,
	0.11457712,
	0.11555660,
	0.11653969,
	0.11752639,
	0.11851670,
	0.11951061,
	0.12050809,
	0.12150916,
	0.12251379,
	0.12352197,
	0.12453371,
	0.12554898,
	0.12656778,
	0.12759009,
	0.12861592,
	0.12964524,
	0.13067806,
	0.13171435,
	0.13275411,
	0.13379734,
	0.13484401,
	0.13589412,
	0.13694767,
	0.13800463,
	0.13906501,
	0.14012878,
	0.14119595,
	0.14226650,
	0.14334041,
	0.14441769,
	0.14549832,
	0.14658229,
	0.14766959,
	0.14876020,
	0.14985413,
	0.15095135,
	0.15205187,
	0.15315566,
	0.15426272,
	0.15537304,
	0.15648660,
	0.15760340,
	0.15872343,
	0.15984667,
	0.16097312,
	0.16210276,
	0.16323558,
	0.16437158,
	0.16551074,
	0.16665305,
	0.16779850,
	0.16894708,
	0.17009878,
	0.17125359,
	0.17241150,
	0.17357249,
	0.17473656,
	0.17590369,
	0.17707388,
	0.17824710,
	0.17942336,
	0.18060264,
	0.18178493,
	0.18297022,
	0.18415849,
	0.18534974,
	0.18654396,
	0.18774113,
	0.18894124,
	0.19014428,
	0.19135024,
	0.19255910,
	0.19377087,
	0.19498552,
	0.19620304,
	0.19742343,
	0.19864666,
	0.19987274,
	0.20110164,
	0.20233336,
	0.20356788,
	0.20480520,
	0.20604530,
	0.20728816,
	0.20853379,
	0.20978216,
	0.21103326,
	0.21228709,
	0.21354363,
	0.21480287,
	0.21606479,
	0.21732939,
	0.21859666,
	0.21986657,
	0.22113912,
	0.22241430,
	0.22369210,
	0.22497250,
	0.22625549,
	0.22754106,
	0.22882919,
	0.23011989,
	0.23141312,
	0.23270888,
	0.23400717,
	0.23530796,
	0.23661124,
	0.23791700,
	0.23922524,
	0.24053593,
	0.24184906,
	0.24316463,
	0.24448262,
	0.24580301,
	0.24712580,
	0.24845097,
	0.24977851,
	0.25110841,
	0.25244066,
	0.25377523,
	0.25511213,
	0.25645133,
	0.25779283,
	0.25913661,
	0.26048266,
	0.26183097,
	0.26318152,
	0.26453430,
	0.26588930,
	0.26724650,
	0.26860590,
	0.26996748,
	0.27133123,
	0.27269713,
	0.27406517,
	0.27543534,
	0.27680763,
	0.27818201,
	0.27955849,
	0.28093705,
	0.28231767,
	0.28370034,
	0.28508505,
	0.28647178,
	0.28786053,
	0.28925127,
	0.29064400,
	0.29203870,
	0.29343536,
	0.29483397,
	0.29623451,
	0.29763697,
	0.29904134,
	0.30044760,
	0.30185574,
	0.30326574,
	0.30467760,
	0.30609130,
	0.30750683,
	0.30892417,
	0.31034331,
	0.31176424,
	0.31318694,
	0.31461141,
	0.31603761,
	0.31746556,
	0.31889522,
	0.32032659,
	0.32175965,
	0.32319439,
	0.32463079,
	0.32606885,
	0.32750855,
	0.32894987,
	0.33039280,
	0.33183733,
	0.33328345,
	0.33473114,
	0.33618038,
	0.33763117,
	0.33908348,
	0.34053732,
	0.34199265,
	0.34344947,
	0.34490777,
	0.34636753,
	0.34782874,
	0.34929138,
	0.35075544,
	0.35222091,
	0.35368777,
	0.35515601,
	0.35662561,
	0.35809656,
	0.35956885,
	0.36104247,
	0.36251739,
	0.36399361,
	0.36547111,
	0.36694988,
	0.36842990,
	0.36991116,
	0.37139364,
	0.37287734,
	0.37436223,
	0.37584831,
	0.37733556,
	0.37882397,
	0.38031351,
	0.38180419,
	0.38329597,
	0.38478886,
	0.38628283,
	0.38777787,
	0.38927398,
	0.39077112,
	0.39226929,
	0.39376848,
	0.39526867,
	0.39676985,
	0.39827199,
	0.39977510,
	0.40127915,
	0.40278413,
	0.40429003,
	0.40579683,
	0.40730451,
	0.40881307,
	0.41032249,
	0.41183276,
	0.41334385,
	0.41485576,
	0.41636847,
	0.41788198,
	0.41939625,
	0.42091128,
	0.42242706,
	0.42394357,
	0.42546080,
	0.42697873,
	0.42849735,
	0.43001664,
	0.43153659,
	0.43305718,
	0.43457841,
	0.43610025,
	0.43762269,
	0.43914572,
	0.44066933,
	0.44219349,
	0.44371820,
	0.44524344,
	0.44676920,
	0.44829545,
	0.44982220,
	0.45134941,
	0.45287709,
	0.45440521,
	0.45593375,
	0.45746272,
	0.45899208,
	0.46052183,
	0.46205195,
	0.46358243,
	0.46511326,
	0.46664441,
	0.46817588,
	0.46970764,
	0.47123969,
	0.47277201,
	0.47430459,
	0.47583741,
	0.47737046,
	0.47890372,
	0.48043718,
	0.48197083,
	0.48350464,
	0.48503861,
	0.48657273,
	0.48810696,
	0.48964131,
	0.49117576,
	0.49271029,
	0.49424489,
	0.49577954,
	0.49731424,
	0.49884895,
	0.50038368,
	0.50191841,
	0.50345312,
	0.50498779,
	0.50652242,
	0.50805698,
	0.50959148,
	0.51112588,
	0.51266017,
	0.51419435,
	0.51572839,
	0.51726228,
	0.51879601,
	0.52032957,
	0.52186293,
	0.52339609,
	0.52492903,
	0.52646173,
	0.52799418,
	0.52952637,
	0.53105828,
	0.53258990,
	0.53412121,
	0.53565220,
	0.53718285,
	0.53871315,
	0.54024309,
	0.54177265,
	0.54330182,
	0.54483057,
	0.54635891,
	0.54788681,
	0.54941425,
	0.55094124,
	0.55246774,
	0.55399374,
	0.55551924,
	0.55704422,
	0.55856866,
	0.56009254,
	0.56161586,
	0.56313860,
	0.56466075,
	0.56618228,
	0.56770320,
	0.56922347,
	0.57074309,
	0.57226205,
	0.57378032,
	0.57529790,
	0.57681477,
	0.57833092,
	0.57984633,
	0.58136098,
	0.58287487,
	0.58438798,
	0.58590029,
	0.58741180,
	0.58892248,
	0.59043232,
	0.59194131,
	0.59344944,
	0.59495668,
	0.59646303,
	0.59796847,
	0.59947299,
	0.60097657,
	0.60247920,
	0.60398087,
	0.60548155,
	0.60698124,
	0.60847992,
	0.60997758,
	0.61147421,
	0.61296978,
	0.61446429,
	0.61595772,
	0.61745006,
	0.61894129,
	0.62043140,
	0.62192038,
	0.62340821,
	0.62489487,
	0.62638036,
	0.62786466,
	0.62934775,
	0.63082963,
	0.63231027,
	0.63378967,
	0.63526780,
	0.63674466,
	0.63822023,
	0.63969450,
	0.64116746,
	0.64263908,
	0.64410936,
	0.64557828,
	0.64704583,
	0.64851200,
	0.64997676,
	0.65144012,
	0.65290204,
	0.65436253,
	0.65582156,
	0.65727912,
	0.65873520,
	0.66018979,
	0.66164287,
	0.66309442,
	0.66454444,
	0.66599290,
	0.66743981,
	0.66888513,
	0.67032887,
	0.67177099,
	0.67321151,
	0.67465038,
	0.67608762,
	0.67752319,
	0.67895709,
	0.68038931,
	0.68181983,
	0.68324863,
	0.68467571,
	0.68610104,
	0.68752463,
	0.68894644,
	0.69036648,
	0.69178472,
	0.69320116,
	0.69461578,
	0.69602856,
	0.69743949,
	0.69884857,
	0.70025577,
	0.70166109,
	0.70306450,
	0.70446600,
	0.70586558,
	0.70726321,
	0.70865890,
	0.71005261,
	0.71144435,
	0.71283410,
	0.71422184,
	0.71560756,
	0.71699125,
	0.71837290,
	0.71975249,
	0.72113001,
	0.72250544,
	0.72387878,
	0.72525001,
	0.72661912,
	0.72798609,
	0.72935092,
	0.73071358,
	0.73207407,
	0.73343237,
	0.73478848,
	0.73614237,
	0.73749404,
	0.73884347,
	0.74019065,
	0.74153556,
	0.74287820,
	0.74421856,
	0.74555661,
	0.74689235,
	0.74822576,
	0.74955683,
	0.75088555,
	0.75221191,
	0.75353589,
	0.75485749,
	0.75617668,
	0.75749346,
	0.75880781,
	0.76011972,
	0.76142919,
	0.76273619,
	0.76404071,
	0.76534275,
	0.76664229,
	0.76793931,
	0.76923381,
	0.77052578,
	0.77181519,
	0.77310205,
	0.77438633,
	0.77566803,
	0.77694713,
	0.77822361,
	0.77949748,
	0.78076872,
	0.78203731,
	0.78330324,
	0.78456650,
	0.78582709,
	0.78708498,
	0.78834016,
	0.78959263,
	0.79084237,
	0.79208937,
	0.79333362,
	0.79457510,
	0.79581381,
	0.79704973,
	0.79828285,
	0.79951316,
	0.80074066,
	0.80196531,
	0.80318712,
	0.80440608,
	0.80562217,
	0.80683538,
	0.80804569,
	0.80925311,
	0.81045761,
	0.81165919,
	0.81285783,
	0.81405352,
	0.81524625,
	0.81643602,
	0.81762280,
	0.81880659,
	0.81998737,
	0.82116514,
	0.82233989,
	0.82351160,
	0.82468026,
	0.82584586,
	0.82700839,
	0.82816784,
	0.82932420,
	0.83047746,
	0.83162760,
	0.83277462,
	0.83391850,
	0.83505924,
	0.83619682,
	0.83733123,
	0.83846246,
	0.83959051,
	0.84071535,
	0.84183699,
	0.84295540,
	0.84407059,
	0.84518253,
	0.84629122,
	0.84739665,
	0.84849880,
	0.84959767,
	0.85069325,
	0.85178552,
	0.85287448,
	0.85396011,
	0.85504241,
	0.85612137,
	0.85719697,
	0.85826920,
	0.85933806,
	0.86040353,
	0.86146561,
	0.86252428,
	0.86357953,
	0.86463136,
	0.86567976,
	0.86672471,
	0.86776620,
	0.86880423,
	0.86983879,
	0.87086986,
	0.87189743,
	0.87292151,
	0.87394206,
	0.87495910,
	0.87597260,
	0.87698256,
	0.87798897,
	0.87899182,
	0.87999110,
	0.88098679,
	0.88197890,
	0.88296741,
	0.88395231,
	0.88493359,
	0.88591125,
	0.88688527,
	0.88785565,
	0.88882237,
	0.88978542,
	0.89074481,
	0.89170051,
	0.89265253,
	0.89360084,
	0.89454545,
	0.89548634,
	0.89642350,
	0.89735693,
	0.89828661,
	0.89921254,
	0.90013471,
	0.90105311,
	0.90196773,
	0.90287857,
	0.90378560,
	0.90468884,
	0.90558826,
	0.90648386,
	0.90737563,
	0.90826356,
	0.90914765,
	0.91002788,
	0.91090425,
	0.91177675,
	0.91264536,
	0.91351009,
	0.91437093,
	0.91522785,
	0.91608087,
	0.91692997,
	0.91777514,
	0.91861637,
	0.91945366,
	0.92028700,
	0.92111637,
	0.92194178,
	0.92276322,
	0.92358067,
	0.92439413,
	0.92520359,
	0.92600905,
	0.92681049,
	0.92760791,
	0.92840130,
	0.92919066,
	0.92997597,
	0.93075723,
	0.93153443,
	0.93230757,
	0.93307663,
	0.93384162,
	0.93460251,
	0.93535932,
	0.93611202,
	0.93686061,
	0.93760508,
	0.93834543,
	0.93908166,
	0.93981374,
	0.94054168,
	0.94126547,
	0.94198511,
	0.94270058,
	0.94341188,
	0.94411900,
	0.94482193,
	0.94552068,
	0.94621523,
	0.94690557,
	0.94759171,
	0.94827362,
	0.94895131,
	0.94962478,
	0.95029401,
	0.95095899,
	0.95161973,
	0.95227621,
	0.95292843,
	0.95357638,
	0.95422006,
	0.95485946,
	0.95549457,
	0.95612539,
	0.95675192,
	0.95737414,
	0.95799205,
	0.95860565,
	0.95921493,
	0.95981988,
	0.96042050,
	0.96101678,
	0.96160871,
	0.96219630,
	0.96277953,
	0.96335841,
	0.96393291,
	0.96450305,
	0.96506881,
	0.96563019,
	0.96618718,
	0.96673978,
	0.96728798,
	0.96783178,
	0.96837117,
	0.96890615,
	0.96943671,
	0.96996285,
	0.97048456,
	0.97100183,
	0.97151467,
	0.97202307,
	0.97252702,
	0.97302652,
	0.97352156,
	0.97401214,
	0.97449826,
	0.97497990,
	0.97545707,
	0.97592976,
	0.97639796,
	0.97686168,
	0.97732090,
	0.97777563,
	0.97822586,
	0.97867158,
	0.97911279,
	0.97954948,
	0.97998166,
	0.98040932,
	0.98083244,
	0.98125104,
	0.98166511,
	0.98207464,
	0.98247962,
	0.98288006,
	0.98327595,
	0.98366729,
	0.98405407,
	0.98443629,
	0.98481394,
	0.98518703,
	0.98555554,
	0.98591949,
	0.98627885,
	0.98663363,
	0.98698383,
	0.98732944,
	0.98767045,
	0.98800688,
	0.98833870,
	0.98866593,
	0.98898855,
	0.98930656,
	0.98961996,
	0.98992875,
	0.99023293,
	0.99053248,
	0.99082742,
	0.99111772,
	0.99140341,
	0.99168446,
	0.99196088,
	0.99223266,
	0.99249981,
	0.99276232,
	0.99302018,
	0.99327340,
	0.99352197,
	0.99376590,
	0.99400517,
	0.99423978,
	0.99446974,
	0.99469504,
	0.99491568,
	0.99513166,
	0.99534297,
	0.99554962,
	0.99575159,
	0.99594890,
	0.99614153,
	0.99632949,
	0.99651277,
	0.99669138,
	0.99686530,
	0.99703455,
	0.99719911,
	0.99735898,
	0.99751417,
	0.99766468,
	0.99781049,
	0.99795162,
	0.99808805,
	0.99821979,
	0.99834683,
	0.99846919,
	0.99858684,
	0.99869980,
	0.99880805,
	0.99891161,
	0.99901047,
	0.99910463,
	0.99919408,
	0.99927883,
	0.99935888,
	0.99943422,
	0.99950486,
	0.99957079,
	0.99963201,
	0.99968853,
	0.99974034,
	0.99978744,
	0.99982983,
	0.99986751,
	0.99990049,
	0.99992875,
	0.99995230,
	0.99997115,
	0.99998528,
	0.99999470,
	0.99999941,
	0.99999941,
	0.99999470,
	0.99998528,
	0.99997115,
	0.99995230,
	0.99992875,
	0.99990049,
	0.99986751,
	0.99982983,
	0.99978744,
	0.99974034,
	0.99968853,
	0.99963201,
	0.99957079,
	0.99950486,
	0.99943422,
	0.99935888,
	0.99927883,
	0.99919408,
	0.99910463,
	0.99901047,
	0.99891161,
	0.99880805,
	0.99869980,
	0.99858684,
	0.99846919,
	0.99834683,
	0.99821979,
	0.99808805,
	0.99795162,
	0.99781049,
	0.99766468,
	0.99751417,
	0.99735898,
	0.99719911,
	0.99703455,
	0.99686530,
	0.99669138,
	0.99651277,
	0.99632949,
	0.99614153,
	0.99594890,
	0.99575159,
	0.99554962,
	0.99534297,
	0.99513166,
	0.99491568,
	0.99469504,
	0.99446974,
	0.99423978,
	0.99400517,
	0.99376590,
	0.99352197,
	0.99327340,
	0.99302018,
	0.99276232,
	0.99249981,
	0.99223266,
	0.99196088,
	0.99168446,
	0.99140341,
	0.99111772,
	0.99082742,
	0.99053248,
	0.99023293,
	0.98992875,
	0.98961996,
	0.98930656,
	0.98898855,
	0.98866593,
	0.98833870,
	0.98800688,
	0.98767045,
	0.98732944,
	0.98698383,
	0.98663363,
	0.98627885,
	0.98591949,
	0.98555554,
	0.98518703,
	0.98481394,
	0.98443629,
	0.98405407,
	0.98366729,
	0.98327595,
	0.98288006,
	0.98247962,
	0.98207464,
	0.98166511,
	0.98125104,
	0.98083244,
	0.98040932,
	0.97998166,
	0.97954948,
	0.97911279,
	0.97867158,
	0.97822586,
	0.97777563,
	0.97732090,
	0.97686168,
	0.97639796,
	0.97592976,
	0.97545707,
	0.97497990,
	0.97449826,
	0.97401214,
	0.97352156,
	0.97302652,
	0.97252702,
	0.97202307,
	0.97151467,
	0.97100183,
	0.97048456,
	0.96996285,
	0.96943671,
	0.96890615,
	0.96837117,
	0.96783178,
	0.96728798,
	0.96673978,
	0.96618718,
	0.96563019,
	0.96506881,
	0.96450305,
	0.96393291,
	0.96335841,
	0.96277953,
	0.96219630,
	0.96160871,
	0.96101678,
	0.96042050,
	0.95981988,
	0.95921493,
	0.95860565,
	0.95799205,
	0.95737414,
	0.95675192,
	0.95612539,
	0.95549457,
	0.95485946,
	0.95422006,
	0.95357638,
	0.95292843,
	0.95227621,
	0.95161973,
	0.95095899,
	0.95029401,
	0.94962478,
	0.94895131,
	0.94827362,
	0.94759171,
	0.94690557,
	0.94621523,
	0.94552068,
	0.94482193,
	0.94411900,
	0.94341188,
	0.94270058,
	0.94198511,
	0.94126547,
	0.94054168,
	0.93981374,
	0.93908166,
	0.93834543,
	0.93760508,
	0.93686061,
	0.93611202,
	0.93535932,
	0.93460251,
	0.93384162,
	0.93307663,
	0.93230757,
	0.93153443,
	0.93075723,
	0.92997597,
	0.92919066,
	0.92840130,
	0.92760791,
	0.92681049,
	0.92600905,
	0.92520359,
	0.92439413,
	0.92358067,
	0.92276322,
	0.92194178,
	0.92111637,
	0.92028700,
	0.91945366,
	0.91861637,
	0.91777514,
	0.91692997,
	0.91608087,
	0.91522785,
	0.91437093,
	0.91351009,
	0.91264536,
	0.91177675,
	0.91090425,
	0.91002788,
	0.90914765,
	0.90826356,
	0.90737563,
	0.90648386,
	0.90558826,
	0.90468884,
	0.90378560,
	0.90287857,
	0.90196773,
	0.90105311,
	0.90013471,
	0.89921254,
	0.89828661,
	0.89735693,
	0.89642350,
	0.89548634,
	0.89454545,
	0.89360084,
	0.89265253,
	0.89170051,
	0.89074481,
	0.88978542,
	0.88882237,
	0.88785565,
	0.88688527,
	0.88591125,
	0.88493359,
	0.88395231,
	0.88296741,
	0.88197890,
	0.88098679,
	0.87999110,
	0.87899182,
	0.87798897,
	0.87698256,
	0.87597260,
	0.87495910,
	0.87394206,
	0.87292151,
	0.87189743,
	0.87086986,
	0.86983879,
	0.86880423,
	0.86776620,
	0.86672471,
	0.86567976,
	0.86463136,
	0.86357953,
	0.86252428,
	0.86146561,
	0.86040353,
	0.85933806,
	0.85826920,
	0.85719697,
	0.85612137,
	0.85504241,
	0.85396011,
	0.85287448,
	0.85178552,
	0.85069325,
	0.84959767,
	0.84849880,
	0.84739665,
	0.84629122,
	0.84518253,
	0.84407059,
	0.84295540,
	0.84183699,
	0.84071535,
	0.83959051,
	0.83846246,
	0.83733123,
	0.83619682,
	0.83505924,
	0.83391850,
	0.83277462,
	0.83162760,
	0.83047746,
	0.82932420,
	0.82816784,
	0.82700839,
	0.82584586,
	0.82468026,
	0.82351160,
	0.82233989,
	0.82116514,
	0.81998737,
	0.81880659,
	0.81762280,
	0.81643602,
	0.81524625,
	0.81405352,
	0.81285783,
	0.81165919,
	0.81045761,
	0.80925311,
	0.80804569,
	0.80683538,
	0.80562217,
	0.80440608,
	0.80318712,
	0.80196531,
	0.80074066,
	0.79951316,
	0.79828285,
	0.79704973,
	0.79581381,
	0.79457510,
	0.79333362,
	0.79208937,
	0.79084237,
	0.78959263,
	0.78834016,
	0.78708498,
	0.78582709,
	0.78456650,
	0.78330324,
	0.78203731,
	0.78076872,
	0.77949748,
	0.77822361,
	0.77694713,
	0.77566803,
	0.77438633,
	0.77310205,
	0.77181519,
	0.77052578,
	0.76923381,
	0.76793931,
	0.76664229,
	0.76534275,
	0.76404071,
	0.76273619,
	0.76142919,
	0.76011972,
	0.75880781,
	0.75749346,
	0.75617668,
	0.75485749,
	0.75353589,
	0.75221191,
	0.75088555,
	0.74955683,
	0.74822576,
	0.74689235,
	0.74555661,
	0.74421856,
	0.74287820,
	0.74153556,
	0.74019065,
	0.73884347,
	0.73749404,
	0.73614237,
	0.73478848,
	0.73343237,
	0.73207407,
	0.73071358,
	0.72935092,
	0.72798609,
	0.72661912,
	0.72525001,
	0.72387878,
	0.72250544,
	0.72113001,
	0.71975249,
	0.71837290,
	0.71699125,
	0.71560756,
	0.71422184,
	0.71283410,
	0.71144435,
	0.71005261,
	0.70865890,
	0.70726321,
	0.70586558,
	0.70446600,
	0.70306450,
	0.70166109,
	0.70025577,
	0.69884857,
	0.69743949,
	0.69602856,
	0.69461578,
	0.69320116,
	0.69178472,
	0.69036648,
	0.68894644,
	0.68752463,
	0.68610104,
	0.68467571,
	0.68324863,
	0.68181983,
	0.68038931,
	0.67895709,
	0.67752319,
	0.67608762,
	0.67465038,
	0.67321151,
	0.67177099,
	0.67032887,
	0.66888513,
	0.66743981,
	0.66599290,
	0.66454444,
	0.66309442,
	0.66164287,
	0.66018979,
	0.65873520,
	0.65727912,
	0.65582156,
	0.65436253,
	0.65290204,
	0.65144012,
	0.64997676,
	0.64851200,
	0.64704583,
	0.64557828,
	0.64410936,
	0.64263908,
	0.64116746,
	0.63969450,
	0.63822023,
	0.63674466,
	0.63526780,
	0.63378967,
	0.63231027,
	0.63082963,
	0.62934775,
	0.62786466,
	0.62638036,
	0.62489487,
	0.62340821,
	0.62192038,
	0.62043140,
	0.61894129,
	0.61745006,
	0.61595772,
	0.61446429,
	0.61296978,
	0.61147421,
	0.60997758,
	0.60847992,
	0.60698124,
	0.60548155,
	0.60398087,
	0.60247920,
	0.60097657,
	0.59947299,
	0.59796847,
	0.59646303,
	0.59495668,
	0.59344944,
	0.59194131,
	0.59043232,
	0.58892248,
	0.58741180,
	0.58590029,
	0.58438798,
	0.58287487,
	0.58136098,
	0.57984633,
	0.57833092,
	0.57681477,
	0.57529790,
	0.57378032,
	0.57226205,
	0.57074309,
	0.56922347,
	0.56770320,
	0.56618228,
	0.56466075,
	0.56313860,
	0.56161586,
	0.56009254,
	0.55856866,
	0.55704422,
	0.55551924,
	0.55399374,
	0.55246774,
	0.55094124,
	0.54941425,
	0.54788681,
	0.54635891,
	0.54483057,
	0.54330182,
	0.54177265,
	0.54024309,
	0.53871315,
	0.53718285,
	0.53565220,
	0.53412121,
	0.53258990,
	0.53105828,
	0.52952637,
	0.52799418,
	0.52646173,
	0.52492903,
	0.52339609,
	0.52186293,
	0.52032957,
	0.51879601,
	0.51726228,
	0.51572839,
	0.51419435,
	0.51266017,
	0.51112588,
	0.50959148,
	0.50805698,
	0.50652242,
	0.50498779,
	0.50345312,
	0.50191841,
	0.50038368,
	0.49884895,
	0.49731424,
	0.49577954,
	0.49424489,
	0.49271029,
	0.49117576,
	0.48964131,
	0.48810696,
	0.48657273,
	0.48503861,
	0.48350464,
	0.48197083,
	0.48043718,
	0.47890372,
	0.47737046,
	0.47583741,
	0.47430459,
	0.47277201,
	0.47123969,
	0.46970764,
	0.46817588,
	0.46664441,
	0.46511326,
	0.46358243,
	0.46205195,
	0.46052183,
	0.45899208,
	0.45746272,
	0.45593375,
	0.45440521,
	0.45287709,
	0.45134941,
	0.44982220,
	0.44829545,
	0.44676920,
	0.44524344,
	0.44371820,
	0.44219349,
	0.44066933,
	0.43914572,
	0.43762269,
	0.43610025,
	0.43457841,
	0.43305718,
	0.43153659,
	0.43001664,
	0.42849735,
	0.42697873,
	0.42546080,
	0.42394357,
	0.42242706,
	0.42091128,
	0.41939625,
	0.41788198,
	0.41636847,
	0.41485576,
	0.41334385,
	0.41183276,
	0.41032249,
	0.40881307,
	0.40730451,
	0.40579683,
	0.40429003,
	0.40278413,
	0.40127915,
	0.39977510,
	0.39827199,
	0.39676985,
	0.39526867,
	0.39376848,
	0.39226929,
	0.39077112,
	0.38927398,
	0.38777787,
	0.38628283,
	0.38478886,
	0.38329597,
	0.38180419,
	0.38031351,
	0.37882397,
	0.37733556,
	0.37584831,
	0.37436223,
	0.37287734,
	0.37139364,
	0.36991116,
	0.36842990,
	0.36694988,
	0.36547111,
	0.36399361,
	0.36251739,
	0.36104247,
	0.35956885,
	0.35809656,
	0.35662561,
	0.35515601,
	0.35368777,
	0.35222091,
	0.35075544,
	0.34929138,
	0.34782874,
	0.34636753,
	0.34490777,
	0.34344947,
	0.34199265,
	0.34053732,
	0.33908348,
	0.33763117,
	0.33618038,
	0.33473114,
	0.33328345,
	0.33183733,
	0.33039280,
	0.32894987,
	0.32750855,
	0.32606885,
	0.32463079,
	0.32319439,
	0.32175965,
	0.32032659,
	0.31889522,
	0.31746556,
	0.31603761,
	0.31461141,
	0.31318694,
	0.31176424,
	0.31034331,
	0.30892417,
	0.30750683,
	0.30609130,
	0.30467760,
	0.30326574,
	0.30185574,
	0.30044760,
	0.29904134,
	0.29763697,
	0.29623451,
	0.29483397,
	0.29343536,
	0.29203870,
	0.29064400,
	0.28925127,
	0.28786053,
	0.28647178,
	0.28508505,
	0.28370034,
	0.28231767,
	0.28093705,
	0.27955849,
	0.27818201,
	0.27680763,
	0.27543534,
	0.27406517,
	0.27269713,
	0.27133123,
	0.26996748,
	0.26860590,
	0.26724650,
	0.26588930,
	0.26453430,
	0.26318152,
	0.26183097,
	0.26048266,
	0.25913661,
	0.25779283,
	0.25645133,
	0.25511213,
	0.25377523,
	0.25244066,
	0.25110841,
	0.24977851,
	0.24845097,
	0.24712580,
	0.24580301,
	0.24448262,
	0.24316463,
	0.24184906,
	0.24053593,
	0.23922524,
	0.23791700,
	0.23661124,
	0.23530796,
	0.23400717,
	0.23270888,
	0.23141312,
	0.23011989,
	0.22882919,
	0.22754106,
	0.22625549,
	0.22497250,
	0.22369210,
	0.22241430,
	0.22113912,
	0.21986657,
	0.21859666,
	0.21732939,
	0.21606479,
	0.21480287,
	0.21354363,
	0.21228709,
	0.21103326,
	0.20978216,
	0.20853379,
	0.20728816,
	0.20604530,
	0.20480520,
	0.20356788,
	0.20233336,
	0.20110164,
	0.19987274,
	0.19864666,
	0.19742343,
	0.19620304,
	0.19498552,
	0.19377087,
	0.19255910,
	0.19135024,
	0.19014428,
	0.18894124,
	0.18774113,
	0.18654396,
	0.18534974,
	0.18415849,
	0.18297022,
	0.18178493,
	0.18060264,
	0.17942336,
	0.17824710,
	0.17707388,
	0.17590369,
	0.17473656,
	0.17357249,
	0.17241150,
	0.17125359,
	0.17009878,
	0.16894708,
	0.16779850,
	0.16665305,
	0.16551074,
	0.16437158,
	0.16323558,
	0.16210276,
	0.16097312,
	0.15984667,
	0.15872343,
	0.15760340,
	0.15648660,
	0.15537304,
	0.15426272,
	0.15315566,
	0.15205187,
	0.15095135,
	0.14985413,
	0.14876020,
	0.14766959,
	0.14658229,
	0.14549832,
	0.14441769,
	0.14334041,
	0.14226650,
	0.14119595,
	0.14012878,
	0.13906501,
	0.13800463,
	0.13694767,
	0.13589412,
	0.13484401,
	0.13379734,
	0.13275411,
	0.13171435,
	0.13067806,
	0.12964524,
	0.12861592,
	0.12759009,
	0.12656778,
	0.12554898,
	0.12453371,
	0.12352197,
	0.12251379,
	0.12150916,
	0.12050809,
	0.11951061,
	0.11851670,
	0.11752639,
	0.11653969,
	0.11555660,
	0.11457712,
	0.11360129,
	0.11262909,
	0.11166054,
	0.11069565,
	0.10973442,
	0.10877688,
	0.10782302,
	0.10687285,
	0.10592639,
	0.10498364,
	0.10404462,
	0.10310932,
	0.10217776,
	0.10124996,
	0.10032590,
	0.09940562,
	0.09848911,
	0.09757638,
	0.09666744,
	0.09576230,
	0.09486097,
	0.09396346,
	0.09306977,
	0.09217992,
	0.09129391,
	0.09041175,
	0.08953345,
	0.08865902,
	0.08778846,
	0.08692179,
	0.08605900,
	0.08520012,
	0.08434515,
	0.08349409,
	0.08264695,
	0.08180375,
	0.08096449,
	0.08012918,
	0.07929782,
	0.07847042,
	0.07764700,
	0.07682756,
	0.07601210,
	0.07520064,
	0.07439318,
	0.07358973,
	0.07279030,
	0.07199489,
	0.07120352,
	0.07041618,
	0.06963289,
	0.06885366,
	0.06807849,
	0.06730739,
	0.06654036,
	0.06577742,
	0.06501857,
	0.06426382,
	0.06351317,
	0.06276664,
	0.06202423,
	0.06128594,
	0.06055178,
	0.05982177,
	0.05909590,
	0.05837419,
	0.05765664,
	0.05694325,
	0.05623404,
	0.05552901,
	0.05482817,
	0.05413152,
	0.05343907,
	0.05275083,
	0.05206681,
	0.05138700,
	0.05071142,
	0.05004008,
	0.04937297,
	0.04871011,
	0.04805150,
	0.04739715,
	0.04674706,
	0.04610125,
	0.04545971,
	0.04482245,
	0.04418948,
	0.04356081,
	0.04293643,
	0.04231636,
	0.04170061,
	0.04108917,
	0.04048206,
	0.03987927,
	0.03928082,
	0.03868671,
	0.03809695,
	0.03751154,
	0.03693048,
	0.03635379,
	0.03578147,
	0.03521352,
	0.03464995,
	0.03409077,
	0.03353597,
	0.03298557,
	0.03243957,
	0.03189797,
	0.03136079,
	0.03082802,
	0.03029967,
	0.02977574,
	0.02925625,
	0.02874119,
	0.02823057,
	0.02772440,
	0.02722267,
	0.02672540,
	0.02623259,
	0.02574424,
	0.02526036,
	0.02478096,
	0.02430603,
	0.02383558,
	0.02336962,
	0.02290815,
	0.02245117,
	0.02199869,
	0.02155072,
	0.02110725,
	0.02066830,
	0.02023386,
	0.01980395,
	0.01937855,
	0.01895769,
	0.01854136,
	0.01812956,
	0.01772230,
	0.01731959,
	0.01692143,
	0.01652781,
	0.01613875,
	0.01575425,
	0.01537432,
	0.01499894,
	0.01462814,
	0.01426191,
	0.01390026,
	0.01354319,
	0.01319070,
	0.01284279,
	0.01249948,
	0.01216076,
	0.01182664,
	0.01149711,
	0.01117219,
	0.01085187,
	0.01053616,
	0.01022507,
	0.00991858,
	0.00961672,
	0.00931947,
	0.00902685,
	0.00873886,
	0.00845549,
	0.00817675,
	0.00790265,
	0.00763318,
	0.00736836,
	0.00710817,
	0.00685263,
	0.00660173,
	0.00635548,
	0.00611389,
	0.00587694,
	0.00564465,
	0.00541702,
	0.00519405,
	0.00497574,
	0.00476210,
	0.00455312,
	0.00434881,
	0.00414917,
	0.00395420,
	0.00376390,
	0.00357828,
	0.00339734,
	0.00322107,
	0.00304949,
	0.00288259,
	0.00272037,
	0.00256283,
	0.00240999,
	0.00226183,
	0.00211836,
	0.00197958,
	0.00184549,
	0.00171610,
	0.00159140,
	0.00147140,
	0.00135610,
	0.00124549,
	0.00113958,
	0.00103837,
	0.00094186,
	0.00085006,
	0.00076295,
	0.00068056,
	0.00060286,
	0.00052987,
	0.00046159,
	0.00039801,
	0.00033914,
	0.00028498,
	0.00023552,
	0.00019077,
	0.00015074,
	0.00011541,
	0.00008479,
	0.00005888,
	0.00003769,
	0.00002120,
	0.00000942,
	0.00000236,
	0.00000000,
};
