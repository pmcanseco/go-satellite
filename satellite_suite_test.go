package satellite

import (
	. "github.com/onsi/ginkgo"
	. "github.com/onsi/gomega"

	"strconv"
	"strings"
	"testing"
)

func TestSatellite(t *testing.T) {
	RegisterFailHandler(Fail)
	RunSpecs(t, "Satellite Suite")
}

var _ = Describe("go-satellite", func() {
	Describe("ParseTLE", func() {
		It("should return correctly parsed values for given ISS#25544", func() {
			sat, err := ParseTLE("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927", "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537", "wgs84")
			Expect(err).ToNot(HaveOccurred())

			Expect(sat.satnum).To(Equal(int64(25544)))
			Expect(sat.epochyr).To(Equal(int64(8)))
			Expect(sat.epochdays).To(Equal(264.51782528))
			Expect(sat.ndot).To(Equal(-2.182e-05))
			Expect(sat.nddot).To(Equal(0.0))
			Expect(sat.bstar).To(Equal(-1.1606e-05))

			Expect(sat.inclo).To(Equal(51.6416))
			Expect(sat.nodeo).To(Equal(247.4627))
			Expect(sat.ecco).To(Equal(0.0006703))
			Expect(sat.argpo).To(Equal(130.536))
			Expect(sat.mo).To(Equal(325.0288))
			Expect(sat.no).To(Equal(15.72125391))
		})

		It("should return correctly parsed values for given NOAA 19#33591", func() {
			sat, err := ParseTLE("1 33591U 09005A   16163.48990228  .00000077  00000-0  66998-4 0  9990", "2 33591  99.0394 120.2160 0013054 232.8317 127.1662 14.12079902378332", "wgs84")
			Expect(err).ToNot(HaveOccurred())

			Expect(sat.satnum).To(Equal(int64(33591)))
			Expect(sat.epochyr).To(Equal(int64(16)))
			Expect(sat.epochdays).To(Equal(163.48990228))
			Expect(sat.ndot).To(Equal(7.7e-7))
			Expect(sat.nddot).To(Equal(0.0))
			Expect(sat.bstar).To(Equal(.66998e-4))

			Expect(sat.inclo).To(Equal(99.0394))
			Expect(sat.nodeo).To(Equal(120.216))
			Expect(sat.ecco).To(Equal(0.0013054))
			Expect(sat.argpo).To(Equal(232.8317))
			Expect(sat.mo).To(Equal(127.1662))
			Expect(sat.no).To(Equal(14.12079902))
		})

		It("should return correctly parsed values for given TITAN 3C#4362", func() {
			sat, err := ParseTLE("1 04632U 70093B   04031.91070959 -.00000084  00000-0  10000-3 0  9955", "2 04632  11.4628 273.1101 1450506 207.6000 143.9350  1.20231981 44145", "wgs84")
			Expect(err).ToNot(HaveOccurred())

			Expect(sat.satnum).To(Equal(int64(4632)))
			Expect(sat.epochyr).To(Equal(int64(4)))
			Expect(sat.epochdays).To(Equal(31.91070959))
			Expect(sat.ndot).To(Equal(-8.4e-7))
			Expect(sat.nddot).To(Equal(0.0))
			Expect(sat.bstar).To(Equal(.1e-3))

			Expect(sat.inclo).To(Equal(11.4628))
			Expect(sat.nodeo).To(Equal(273.1101))
			Expect(sat.ecco).To(Equal(0.1450506))
			Expect(sat.argpo).To(Equal(207.6))
			Expect(sat.mo).To(Equal(143.935))
			Expect(sat.no).To(Equal(1.20231981))
		})
	})

	Describe("Propagate", func() {
		testCases := [8]PropagationTestCase{
			// PropagationTestCase{
			// 	line1: "",
			// 	line2: "",
			// 	grav:  "wgs72",
			// 	testData: ``,
			// },
			PropagationTestCase{
				line1: "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753",
				line2: "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667",
				grav:  "wgs72",
				testData: `0.00000000 7022.46529266 -1400.08296755 0.03995155 1.893841015 6.405893759 4.534807250
360.00000000 -7154.03120202 -3783.17682504 -3536.19412294 4.741887409 -4.151817765 -2.093935425 2000 6 28 0:50:19.733571
720.00000000 -7134.59340119 6531.68641334 3260.27186483 -4.113793027 -2.911922039 -2.557327851 2000 6 28 6:50:19.733571
1080.00000000 5568.53901181 4492.06992591 3863.87641983 -4.209106476 5.159719888 2.744852980 2000 6 28 12:50:19.733571
1440.00000000 -938.55923943 -6268.18748831 -4294.02924751 7.536105209 -0.427127707 0.989878080 2000 6 28 18:50:19.733571
1800.00000000 -9680.56121728 2802.47771354 124.10688038 -0.905874102 -4.659467970 -3.227347517 2000 6 29 0:50:19.733571
2160.00000000 190.19796988 7746.96653614 5110.00675412 -6.112325142 1.527008184 -0.139152358 2000 6 29 6:50:19.733571
2520.00000000 5579.55640116 -3995.61396789 -1518.82108966 4.767927483 5.123185301 4.276837355 2000 6 29 12:50:19.733571
2880.00000000 -8650.73082219 -1914.93811525 -3007.03603443 3.067165127 -4.828384068 -2.515322836 2000 6 29 18:50:19.733571
3240.00000000 -5429.79204164 7574.36493792 3747.39305236 -4.999442110 -1.800561422 -2.229392830 2000 6 30 0:50:19.733571
3600.00000000 6759.04583722 2001.58198220 2783.55192533 -2.180993947 6.402085603 3.644723952 2000 6 30 6:50:19.733571
3960.00000000 -3791.44531559 -5712.95617894 -4533.48630714 6.668817493 -2.516382327 -0.082384354 2000 6 30 12:50:19.733571
4320.00000000 -9060.47373569 4658.70952502 813.68673153 -2.232832783 -4.110453490 -3.157345433 2000 6 30 18:50:19.733571 `,
			},
			PropagationTestCase{
				line1: "1 04632U 70093B   04031.91070959 -.00000084  00000-0  10000-3 0  9955",
				line2: "2 04632  11.4628 273.1101 1450506 207.6000 143.9350  1.20231981 44145",
				grav:  "wgs72",
				testData: `0.00000000 2334.11450085 -41920.44035349 -0.03867437 2.826321032 -0.065091664 0.570936053
-5184.00000000 -29020.02587128 13819.84419063 -5713.33679183 -1.768068390 -3.235371192 -0.395206135 2004 1 28 7:27:25.308584
-5064.00000000 -32982.56870101 -11125.54996609 -6803.28472771 0.617446996 -3.379240041 0.085954707 2004 1 28 9:27:25.308597
-4944.00000000 -22097.68730513 -31583.13829284 -4836.34329328 2.230597499 -2.166594667 0.426443070 2004 1 28 11:27:25.308611
-4896.00000000 -15129.94694545 -36907.74526221 -3487.56256701 2.581167187 -1.524204737 0.504805763 2004 1 28 12:15:25.308600`,
			},
			PropagationTestCase{
				line1: "1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985",
				line2: "2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774",
				grav:  "wgs72",
				testData: `0.00000000 3988.31022699 5498.96657235 0.90055879 -3.290032738 2.357652820 6.496623475
120.00000000 -3935.69800083 409.10980837 5471.33577327 -3.374784183 -6.635211043 -1.942056221 2006 6 25 21:46:43.980124
240.00000000 -1675.12766915 -5683.30432352 -3286.21510937 5.282496925 1.508674259 -5.354872978 2006 6 25 23:46:43.980097
360.00000000 4993.62642836 2890.54969900 -3600.40145627 0.347333429 5.707031557 5.070699638 2006 6 26 1:46:43.980111
480.00000000 -1115.07959514 4015.11691491 5326.99727718 -5.524279443 -4.765738774 2.402255961 2006 6 26 3:46:43.980124
600.00000000 -4329.10008198 -5176.70287935 409.65313857 2.858408303 -2.933091792 -6.509690397 2006 6 26 5:46:43.980097
720.00000000 3692.60030028 -976.24265255 -5623.36447493 3.897257243 6.415554948 1.429112190 2006 6 26 7:46:43.980111
840.00000000 2301.83510037 5723.92394553 2814.61514580 -5.110924966 -0.764510559 5.662120145 2006 6 26 9:46:43.980124
960.00000000 -4990.91637950 -2303.42547880 3920.86335598 -0.993439372 -5.967458360 -4.759110856 2006 6 26 11:46:43.980097
1080.00000000 642.27769977 -4332.89821901 -5183.31523910 5.720542579 4.216573838 -2.846576139 2006 6 26 13:46:43.980111
1200.00000000 4719.78335752 4798.06938996 -943.58851062 -2.294860662 3.492499389 6.408334723 2006 6 26 15:46:43.980124
1320.00000000 -3299.16993602 1576.83168320 5678.67840638 -4.460347074 -6.202025196 -0.885874586 2006 6 26 17:46:43.980097
1440.00000000 -2777.14682335 -5663.16031708 -2462.54889123 4.915493146 0.123328992 -5.896495091 2006 6 26 19:46:43.980111
1560.00000000 4992.31573893 1716.62356770 -4287.86065581 1.640717189 6.071570434 4.338797931 2006 6 26 21:46:43.980124
1680.00000000 -8.22384755 4662.21521668 4905.66411857 -5.891011274 -3.593173872 3.365100460 2006 6 26 23:46:43.980097
1800.00000000 -4966.20137963 -4379.59155037 1349.33347502 1.763172581 -3.981456387 -6.343279443 2006 6 27 1:46:43.980111
1920.00000000 2954.49390331 -2080.65984650 -5754.75038057 4.895893306 5.858184322 0.375474825 2006 6 27 3:46:43.980124
2040.00000000 3363.28794321 5559.55841180 1956.05542266 -4.587378863 0.591943403 6.107838605 2006 6 27 5:46:43.980097
2160.00000000 -4856.66780070 -1107.03450192 4557.21258241 -2.304158557 -6.186437070 -3.956549542 2006 6 27 7:46:43.980111
2280.00000000 -497.84480071 -4863.46005312 -4700.81211217 5.960065407 2.996683369 -3.767123329 2006 6 27 9:46:43.980124
2400.00000000 5241.61936096 3910.75960683 -1857.93473952 -1.124834806 4.406213160 6.148161299 2006 6 27 11:46:43.980097 
2520.00000000 -2451.38045953 2610.60463261 5729.79022069 -5.366560525 -5.500855666 0.187958716 2006 6 27 13:46:43.980111
2640.00000000 -3791.87520638 -5378.82851382 -1575.82737930 4.266273592 -1.199162551 -6.276154080 2006 6 27 15:46:43.980124
2760.00000000 4730.53958356 524.05006433 -4857.29369725 2.918056288 6.135412849 3.495115636 2006 6 27 17:46:43.980097
2880.00000000 1159.27802897 5056.60175495 4353.49418579 -5.968060341 -2.314790406 4.230722669 2006 6 27 19:46:43.980111`,
			},
			PropagationTestCase{
				line1: "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    8",
				line2: "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  105",
				grav:  "wgs72",
				testData: `0.00000000 2328.96975262 -5995.22051338 1719.97297192 2.912073281 -0.983417956 -7.090816210
120.00000000 1020.69234558 2286.56260634 -6191.55565927 -3.746543902 6.467532721 1.827985678
240.00000000 -3226.54349155 3503.70977525 4532.80979343 1.000992116 -5.788042888 5.162585826
360.00000000 2456.10706533 -6071.93855503 1222.89768554 2.679390040 -0.448290811 -7.228792155
480.00000000 787.16457349 2719.91800946 -6043.86662024 -3.759883839 6.277439314 2.397897864
600.00000000 -3110.97648029 3121.73026235 4878.15217035 1.244916056 -6.124880425 4.700576353
720.00000000 2567.56229695 -6112.50383922 713.96374435 2.440245751 0.098109002 -7.319959258
840.00000000 556.05661780 3144.52288201 -5855.34636178 -3.754660143 6.044752775 2.957941672
960.00000000 -2982.47940539 2712.61663711 5192.32330472 1.475566773 -6.427737014 4.202420227
1080.00000000 2663.08964352 -6115.48290885 196.40072866 2.196121564 0.652415093 -7.362824152
1200.00000000 328.54999674 3557.09490552 -5626.21427211 -3.731193288 5.769341172 3.504058731
1320.00000000 -2842.06876757 2278.42343492 5472.33437150 1.691852635 -6.693216335 3.671022712
1440.00000000 2742.55398832 -6079.67009123 -326.39012649 1.948497651 1.211072678 -7.356193131`,
			},
			PropagationTestCase{
				line1: "1 04632U 70093B   04031.91070959 -.00000084  00000-0  10000-3 0  9955",
				line2: "2 04632  11.4628 273.1101 1450506 207.6000 143.9350  1.20231981 44145",
				grav:  "wgs72",
				testData: `0.00000000 2334.11450085 -41920.44035349 -0.03867437 2.826321032 -0.065091664 0.570936053
-5184.00000000 -29020.02587128 13819.84419063 -5713.33679183 -1.768068390 -3.235371192 -0.395206135
-5064.00000000 -32982.56870101 -11125.54996609 -6803.28472771 0.617446996 -3.379240041 0.085954707
-4944.00000000 -22097.68730513 -31583.13829284 -4836.34329328 2.230597499 -2.166594667 0.426443070
-4896.00000000 -15129.94694545 -36907.74526221 -3487.56256701 2.581167187 -1.524204737 0.504805763`,
			},
			PropagationTestCase{
				line1: "1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600",
				line2: "2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119",
				grav:  "wgs72",
				testData: `0.00000000 7534.10987189 41266.39266843 -0.10801028 -3.027168008 0.558848996 0.207982755
120.00000000 -14289.19940414 39469.05530051 1428.62838591 -2.893205245 -1.045447840 0.179634249
240.00000000 -32222.92014955 26916.25425799 2468.59996594 -1.973007929 -2.359335071 0.102539376
360.00000000 -41413.95109398 7055.51656639 2838.90906671 -0.521665080 -3.029172207 -0.002066843
480.00000000 -39402.72251896 -14716.42475223 2441.32678358 1.066928187 -2.878714619 -0.105865729
600.00000000 -26751.08889828 -32515.13982431 1384.38865570 2.366228869 -1.951032799 -0.181018498
720.00000000 -6874.77975542 -41530.38329422 -46.60245459 3.027415087 -0.494671177 -0.207337260
840.00000000 14859.52039042 -39302.58907247 -1465.02482524 2.869609883 1.100123969 -0.177514425
960.00000000 32553.14863770 -26398.88401807 -2485.45866002 1.930064459 2.401574539 -0.099250520
1080.00000000 41365.67576837 -6298.09965811 -2828.05254033 0.459741276 3.051680214 0.006431872
1200.00000000 38858.83295070 15523.39314924 -2396.86850752 -1.140211488 2.867567143 0.110637217
1320.00000000 25701.46068162 33089.42617648 -1308.68556638 -2.428713821 1.897381431 0.184605907
1440.00000000 5501.08137100 41590.27784405 138.32522930 -3.050691874 0.409203052 0.207958133`,
			},
			PropagationTestCase{
				line1: "1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985",
				line2: "2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774",
				grav:  "wgs72",
				testData: `0.00000000 3988.31022699 5498.96657235 0.90055879 -3.290032738 2.357652820 6.496623475
120.00000000 -3935.69800083 409.10980837 5471.33577327 -3.374784183 -6.635211043 -1.942056221
240.00000000 -1675.12766915 -5683.30432352 -3286.21510937 5.282496925 1.508674259 -5.354872978
360.00000000 4993.62642836 2890.54969900 -3600.40145627 0.347333429 5.707031557 5.070699638
480.00000000 -1115.07959514 4015.11691491 5326.99727718 -5.524279443 -4.765738774 2.402255961
600.00000000 -4329.10008198 -5176.70287935 409.65313857 2.858408303 -2.933091792 -6.509690397
720.00000000 3692.60030028 -976.24265255 -5623.36447493 3.897257243 6.415554948 1.429112190
840.00000000 2301.83510037 5723.92394553 2814.61514580 -5.110924966 -0.764510559 5.662120145
960.00000000 -4990.91637950 -2303.42547880 3920.86335598 -0.993439372 -5.967458360 -4.759110856
1080.00000000 642.27769977 -4332.89821901 -5183.31523910 5.720542579 4.216573838 -2.846576139
1200.00000000 4719.78335752 4798.06938996 -943.58851062 -2.294860662 3.492499389 6.408334723
1320.00000000 -3299.16993602 1576.83168320 5678.67840638 -4.460347074 -6.202025196 -0.885874586
1440.00000000 -2777.14682335 -5663.16031708 -2462.54889123 4.915493146 0.123328992 -5.896495091
1560.00000000 4992.31573893 1716.62356770 -4287.86065581 1.640717189 6.071570434 4.338797931
1680.00000000 -8.22384755 4662.21521668 4905.66411857 -5.891011274 -3.593173872 3.365100460
1800.00000000 -4966.20137963 -4379.59155037 1349.33347502 1.763172581 -3.981456387 -6.343279443
1920.00000000 2954.49390331 -2080.65984650 -5754.75038057 4.895893306 5.858184322 0.375474825
2040.00000000 3363.28794321 5559.55841180 1956.05542266 -4.587378863 0.591943403 6.107838605
2160.00000000 -4856.66780070 -1107.03450192 4557.21258241 -2.304158557 -6.186437070 -3.956549542
2280.00000000 -497.84480071 -4863.46005312 -4700.81211217 5.960065407 2.996683369 -3.767123329
2400.00000000 5241.61936096 3910.75960683 -1857.93473952 -1.124834806 4.406213160 6.148161299
2520.00000000 -2451.38045953 2610.60463261 5729.79022069 -5.366560525 -5.500855666 0.187958716
2640.00000000 -3791.87520638 -5378.82851382 -1575.82737930 4.266273592 -1.199162551 -6.276154080
2760.00000000 4730.53958356 524.05006433 -4857.29369725 2.918056288 6.135412849 3.495115636
2880.00000000 1159.27802897 5056.60175495 4353.49418579 -5.968060341 -2.314790406 4.230722669`,
			},
			PropagationTestCase{
				line1: "1 23599U 95029B   06171.76535463  .00085586  12891-6  12956-2 0  2905",
				line2: "2 23599   6.9327   0.2849 5782022 274.4436  25.2425  4.47796565123555",
				grav:  "wgs72",
				testData: `0.00000000 9892.63794341 35.76144969 -1.08228838 3.556643237 6.456009375 0.783610890       
20.00000000 11931.95642997 7340.74973750 886.46365987 0.308329116 5.532328972 0.672887281
40.00000000 11321.71039205 13222.84749156 1602.40119049 -1.151973982 4.285810871 0.521919425
60.00000000 9438.29395675 17688.05450261 2146.59293402 -1.907904054 3.179955046 0.387692479
80.00000000 6872.08634639 20910.11016811 2539.79945034 -2.323995367 2.207398462 0.269506121
100.00000000 3933.37509798 23024.07662542 2798.25966746 -2.542860616 1.327134966 0.162450076
120.00000000 816.64091546 24118.98675475 2932.69459428 -2.626838010 0.504502763 0.062344306
140.00000000 -2334.41705804 24246.86096326 2949.36448841 -2.602259646 -0.288058266 -0.034145135
160.00000000 -5394.31798039 23429.42716149 2850.86832586 -2.474434068 -1.074055982 -0.129868366
180.00000000 -8233.35130237 21661.24480883 2636.51456118 -2.230845533 -1.875742344 -0.227528603
200.00000000 -10693.96497348 18909.88168891 2302.33707548 -1.835912433 -2.716169865 -0.329931880
220.00000000 -12553.89669904 15114.63990716 1840.93573231 -1.212478879 -3.619036996 -0.439970633
240.00000000 -13450.20591864 10190.57904289 1241.95958736 -0.189082511 -4.596701971 -0.559173899
260.00000000 -12686.60437121 4079.31106161 498.27078614 1.664498211 -5.559889865 -0.676747779
280.00000000 -8672.55867753 -2827.56823315 -342.59644716 5.515079852 -5.551222962 -0.676360044
300.00000000 1153.31498060 -6411.98692060 -779.87288941 9.689818102 1.388598425 0.167868798
320.00000000 9542.79201056 -533.71253081 -65.73165428 3.926947087 6.459583539 0.785686755
340.00000000 11868.80960100 6861.59590848 833.72780602 0.452957852 5.632811328 0.685262323
360.00000000 11376.23941678 12858.97121366 1563.40660172 -1.087665695 4.374693347 0.532207051
380.00000000 9547.70300782 17421.48570758 2118.56907515 -1.876540262 3.253891728 0.395810243
400.00000000 7008.51470263 20725.47471227 2520.56064289 -2.308703599 2.270724438 0.276138613
420.00000000 4082.28135104 22911.04184601 2786.37568309 -2.536665546 1.383670232 0.168153407
440.00000000 969.17978149 24071.23673676 2927.31326579 -2.626695115 0.557172428 0.067536854
460.00000000 -2184.71515444 24261.21671601 2950.08142825 -2.607072866 -0.236887607 -0.029125215
480.00000000 -5253.42223370 23505.37595671 2857.66120738 -2.484424544 -1.022255436 -0.124714444
500.00000000 -8108.27961017 21800.81688388 2649.72981961 -2.247597251 -1.821159176 -0.221925624
520.00000000 -10594.77795556 19117.80779221 2322.72136979 -1.863118484 -2.656426668 -0.323521502
540.00000000 -12497.32045995 15398.64085906 1869.69983897 -1.258130763 -3.551583368 -0.432338888
560.00000000 -13467.92475245 10560.90147785 1280.78399181 -0.271870523 -4.520514224 -0.550016092
580.00000000 -12848.18843590 4541.21901842 548.53826427 1.494157156 -5.489585384 -0.667472039
600.00000000 -9152.70552728 -2344.24950144 -287.98121970 5.127921095 -5.650383025 -0.685989008
620.00000000 280.38490909 -6500.10264018 -790.36092984 9.779619614 0.581815811 0.074171345
640.00000000 9166.25784315 -1093.12552651 -129.49428887 4.316668714 6.438636494 0.785116609
660.00000000 11794.48942915 6382.21138354 780.88439015 0.604412453 5.731729369 0.697574333
680.00000000 11424.30138324 12494.26088864 1524.33165488 -1.021328075 4.463448968 0.542532698
700.00000000 9652.09867350 17153.84762075 2090.48038336 -1.844516637 3.327522235 0.403915232
720.00000000 7140.41945884 20539.25485336 2501.21469368 -2.293173684 2.333507912 0.282716311`,
			},
		}

		for _, testCase := range testCases {
			propagationTest(testCase)
		}
	})
})

type PropagationTestCase struct {
	line1, line2, testData string
	grav                   Gravity
}

func propagationTest(testCase PropagationTestCase) {
	satrec, err := TLEToSat(testCase.line1, testCase.line2, testCase.grav)
	if err != nil {
		panic(err)
	}

	lines := strings.Split(testCase.testData, "\n")

	for _, line := range lines {
		Context("Satnum "+strconv.FormatInt(satrec.satnum, 10), func() {
			theoData := strings.Split(line, " ")

			theoPos := Vector3{X: parseFloat(theoData[1]), Y: parseFloat(theoData[2]), Z: parseFloat(theoData[3])}
			theoVel := Vector3{X: parseFloat(theoData[4]), Y: parseFloat(theoData[5]), Z: parseFloat(theoData[6])}

			expPos, expVel := sgp4(satrec, parseFloat(theoData[0]))

			It("Should produce accurate results for time "+theoData[0], func() {
				Expect(expPos.X).To(BeNumerically("~", theoPos.X, 0.0001))
				Expect(expPos.Y).To(BeNumerically("~", theoPos.Y, 0.0001))
				Expect(expPos.Z).To(BeNumerically("~", theoPos.Z, 0.0001))

				Expect(expVel.X).To(BeNumerically("~", theoVel.X, 0.0001))
				Expect(expVel.Y).To(BeNumerically("~", theoVel.Y, 0.0001))
				Expect(expVel.Z).To(BeNumerically("~", theoVel.Z, 0.0001))
			})
		})
	}
}
