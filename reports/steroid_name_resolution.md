# Steroid name resolution: report only, nothing assigned

Produced by `tools/resolve_steroid_names.py` against the 87 distinct compound names in
`data/seed/steroid_jasms2022.csv`, read-only. **No assignment here has entered the corpus.**
Nothing in this file is covered by `corpus_sha256`, and running the script cannot change
any answer the API gives.

A structure assignment is a curation act needing a provenance trail. A name that resolved
through a public service on one day is not one.

**87 names, 6 resolved clean, 67 ambiguous, 14 failed - manual chemist review needed for 81.**

2 of the 6 clean resolutions are additionally marked **NEEDS CONFIRMATION**: they resolved from their systematic name AND share connectivity with another name in this corpus. See below.

| name | InChIKey or FAILED | source | ambiguity |
| --- | --- | --- | --- |
| `1,3,5(10)-estratiene-3,16,17-triol-3-sulfate` | ZORQMBLUMWNJEQ-ZXXIGWHRSA-N | commercial name only | name-only |
| `1,3,5(10)-estratiene-3,17β-diol 17-sulfate` | JSUDNGPWAXYETN-ZBRFXRBCSA-N | commercial name only | name-only |
| `1,3,5(10)-estratiene-3,17β-diol 3-glucuronide` | MUOHJTRCBBDUOW-QXYWQCSFSA-N | commercial name only | name-only |
| `1,3,5(10)-estratiene-3,17β-diol diglucuronide` | FAILED | - |  |
| `1,3,5(10)-estratiene-3-ol-11,17-dione` | WKUQOYYYFVBROI-UHFFFAOYSA-N | commercial name only | name-only |
| `1,3,5(10)-estratiene-3-ol-17-one 3-sulfate` | JKKFKPJIXZFSSB-CBZIJGRNSA-N, JKKFKPJIXZFSSB-UHFFFAOYSA-N | commercial name, resolvers disagree | name-only; stereo |
| `1,3,5(10)-estratriene-3,16α-17β-triol-3-glucuronide` | UZKIAJMSMKLBQE-JRSYHJKYSA-N, UZKIAJMSMKLBQE-WTSDUJKYSA-N | commercial name, resolvers disagree | name-only; stereo |
| `1,3,5(10)-estratriene-3,16α-diol-17-one` | WPOCIZJTELRQMF-LVQHMEKZSA-N, WPOCIZJTELRQMF-QFXBJFAPSA-N, WPOCIZJTELRQMF-UHFFFAOYSA-N | commercial name, resolvers disagree | name-only; stereo |
| `1,3,5(10)-estratriene-3,17α-diol-3-sulfate` | QZIGLSSUDXBTLJ-SFFUCWETSA-N | CIR | **NEEDS CONFIRMATION** - stereo-pair: 2 names share connectivity QZIGLSSUDXBTLJ |
| `1,3,5(10)-estratriene-3,17β-diol-17-glucuronide` | FAILED | - |  |
| `1,3,5(10)-estratriene-3,17β-diol-3-sulfate` | QZIGLSSUDXBTLJ-ZBRFXRBCSA-N | CIR | **NEEDS CONFIRMATION** - stereo-pair: 2 names share connectivity QZIGLSSUDXBTLJ |
| `1,4-androstadiene-17α-methyl-17β-ol-3-one` | GCKMFJBGXUYNAG-HLXURNFRSA-N, XWALNWXLMVGSFR-HLXURNFRSA-N | commercial name, resolvers disagree | name-only; stereo-pair: 2 names share connectivity GCKMFJBGXUYNAG |
| `1,4-androstadiene-17β-ol-3-one` | RSIHSRDYCUFFLA-DYKIIFRCSA-N | commercial name only | name-only |
| `1,4-androstadiene-17β-ol-3-one glucuronide` | WTOYIEKXMIILRQ-UHFFFAOYSA-N | commercial name only | name-only |
| `1,4-androstadiene-17β-ol-3-one propionate` | ULJOJMSGJSWPSE-BLQWBTBKSA-N | commercial name only | name-only |
| `1,4-androstadiene-17β-ol-3-one sulfate` | GSLOLTKGVZFNKZ-DYKIIFRCSA-N, GSLOLTKGVZFNKZ-KZYORJDKSA-N, GSLOLTKGVZFNKZ-UHFFFAOYSA-N | commercial name, resolvers disagree | name-only; stereo |
| `1,4-androstadiene-17β-ol-3-one undecylenate*` | FAILED | - |  |
| `1,4-androstadiene-3,17-dione` | LUJVUUWNAPIQQI-QAGGRKNESA-N | PubChem+CIR |  |
| `1,4-androstadiene-4-chloro-17α-methyl-17β-ol-3-one` | AGUNEISBPXQOPA-XMUHMHRVSA-N | commercial name only | name-only |
| `1,4-pregnadiene-11β,16α,17,21-tetrol-3,20-dione 16,17 acetonide` | WBGKWQHBNHJJPZ-LECWWXJVSA-N, WBGKWQHBNHJJPZ-YWZQBGSISA-N | commercial name, resolvers disagree | name-only; stereo |
| `1,4-pregnadiene-11β,17,21-triol-3,20-dione` | OIGNJSKKLXVSLS-VWUMJDOOSA-N | commercial name only | name-only |
| `1,4-pregnadiene-11β,17,21-triol-3,20-dione 21-hemisuccinate` | APGDTXUMTIZLCJ-CGVGKPPMSA-N | commercial name only | name-only |
| `1,4-pregnadiene-16α-methyl-6α,9-difluoro-11β,21-diol-3,20-dione 21-pivalate` | UWGRWFCLGQWKPR-GSTUPEFVSA-N | commercial name only | name-only |
| `1,4-pregnadiene-17,21-diol-3,11,20-trione` | XOFYZVNMUHMLCC-ZPOLXVRWSA-N | commercial name only | name-only |
| `1,4-pregnadiene-6-fluoro-11β,16α,17,21-tetrol-3,20-dione acetonide` | XSFJVAJPIHIPKU-XWCQMRHXSA-N | commercial name only | name-only |
| `1,4-pregnadiene-6α,9α-difluoro-16α-methyl-11β,17,21-triol-3,20-dione` | WXURHACBFYSXBI-GQKYHHCASA-N | commercial name only | name-only |
| `1,4-pregnadiene-6α-methyl-11β,17α,21-triol-3,20-one` | VHRSUDSXCMQTMA-PJHHCJLFSA-N, VHRSUDSXCMQTMA-UHFFFAOYSA-N | commercial name, resolvers disagree | name-only; stereo |
| `1,4-pregnadiene-9α-fluoro-11β,16α,17,21-tetrol-3,20-dione` | GFNANZIMVAIWHM-OBYCQNJPSA-N | commercial name only | name-only |
| `1,4-pregnadiene-9α-fluoro-11β,16α,17,21-tetrol-3,20-dione 16,21-diacetate` | XGMPVBXKDAHORN-RBWIMXSLSA-N | commercial name only | name-only |
| `1,4-pregnadiene-9α-fluoro-16α-methyl-11β,17,21-triol-3,20-dione` | UREBDLICKHMUKA-CXSFZGCWSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity UREBDLICKHMUKA |
| `1,4-pregnadiene-9α-fluoro-16β-methyl-11β,17,21-triol-3,20-dione` | UREBDLICKHMUKA-DVTGEIKXSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity UREBDLICKHMUKA |
| `17-α-2,4-pregnadiene-20-yno(2,3-d)isoxazol-17-ol` | POZRVZJJTULAOH-LHZXLZLDSA-N | commercial name only | name-only |
| `4,6-pregnadiene-6-methyl-16-methylene-17-ol-3,20-dione acetate` | UDKABVSQKJNZBH-DWNQPYOZSA-N | commercial name only | name-only |
| `4,9,11-estratiene-17α-ol-3-one` | FAILED | - | commercial name names two things |
| `4,9,11-estratiene-17β-ol-3-one` | MEHHPFQKXOUFFV-OWSLCNJRSA-N | commercial name only | name-only |
| `4,9,11-estratrien-13-ethyl-17-ethinyl-17β-ol-3-one` | BJJXHLWLUDYTGC-ANULTFPQSA-N | commercial name only | name-only |
| `4,9,11-estratriene-17β-ol-3-one acetate` | CMRJPMODSSEAPL-FYQPLNBISA-N | commercial name only | name-only |
| `4-androstene-14α,17β-diol-3-one` | FAILED | - |  |
| `4-androstene-16α-17β-diol-3-one` | YMCWOAZGWMZGQT-FPNLOETNSA-N | commercial name only | name-only |
| `4-androstene-17α-methyl-17β-ol-3-one` | GCKMFJBGXUYNAG-HLXURNFRSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity GCKMFJBGXUYNAG |
| `4-androstene-17α-ol-3-one` | MUMGGOZAMZWBJJ-KZYORJDKSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity MUMGGOZAMZWBJJ |
| `4-androstene-17α-ol-3-one glucuronide` | NIKZPECGCSUSBV-HMAFJQTKSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity NIKZPECGCSUSBV |
| `4-androstene-17α-ol-3-one sulfate` | WAQBISPOEAOCOG-FZPSTPAASA-N, WAQBISPOEAOCOG-KZYORJDKSA-N, WAQBISPOEAOCOG-UHFFFAOYSA-N | commercial name, resolvers disagree | name-only; stereo; stereo-pair: 2 names share connectivity WAQBISPOEAOCOG |
| `4-androstene-17β-ol-3-one` | MUMGGOZAMZWBJJ-DYKIIFRCSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity MUMGGOZAMZWBJJ |
| `4-androstene-17β-ol-3-one 17-glucuronide` | NIKZPECGCSUSBV-HMAFJQTKSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity NIKZPECGCSUSBV |
| `4-androstene-17β-ol-3-one acetate` | DJPZSBANTAQNFN-PXQJOHHUSA-N | commercial name only | name-only |
| `4-androstene-17β-ol-3-one benzoate` | RZJSCADWIWNGKI-IXKNJLPQSA-N | commercial name only | name-only |
| `4-androstene-17β-ol-3-one cypionate` | HPFVBGJFAYZEBE-ZLQWOROUSA-N | commercial name only | name-only |
| `4-androstene-17β-ol-3-one enanthate` | VOCBWIIFXDYGNZ-IXKNJLPQSA-N | commercial name only | name-only |
| `4-androstene-17β-ol-3-one isocaproate` | PPYHLSBUTAPNGT-BKWLFHPQSA-N, PPYHLSBUTAPNGT-UHFFFAOYSA-N | commercial name, resolvers disagree | name-only; stereo |
| `4-androstene-17β-ol-3-one sulfate` | WAQBISPOEAOCOG-DYKIIFRCSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity WAQBISPOEAOCOG |
| `4-androstene-3,17-dione` | AEMFNILZOJDQLW-QAGGRKNESA-N | PubChem |  |
| `4-androstene-4-chloro-17α-methyl-17β-ol-3-one` | SOMOGWLYTLQJGT-UHFFFAOYSA-N | commercial name only | name-only |
| `4-androstene-4-chloro-17β-ol-3-one` | FAILED | - | commercial name names two things |
| `4-androstene-4-chloro-17β-ol-3-one acetate` | FAILED | - | commercial name names two things |
| `4-androstene-9α-fluoro-11β-ol-17α-methyl-17β-ol-3-one` | YLRFCQOZQXIBAB-RBZZARIASA-N | commercial name only | name-only |
| `4-estrene-17α-ethinyl-17β-ol-3-one` | FAILED | - | commercial name names two things |
| `4-estrene-17α-ethinyl-17β-ol-3-one acetate` | IMONTRJLAWHYGT-ZCPXKWAGSA-N | commercial name only | name-only |
| `4-estrene-17α-ol-3-one` | FAILED | - |  |
| `4-estrene-17β-ol-3-one` | FAILED | - | commercial name names two things |
| `4-estrene-17β-ol-3-one D-glucuronide` | ISBYSZZUCBXGIH-UHFFFAOYSA-N | commercial name only | name-only |
| `4-estrene-17β-ol-3-one benzoate` | FAILED | - |  |
| `4-estrene-17β-ol-3-one sulfate` | SKZMVWBZTQNCKW-IZPLOLCNSA-N | commercial name only | name-only |
| `4-estrene-4-fluoro-17β-ol-3-one` | KIUHHFFZDPMQQM-YGRHGMIBSA-N | commercial name only | name-only |
| `4-pregnene-11β,17,18,21-tetrol-3,20-dione` | HESFZGWRDUVOMS-UKSDXMLSSA-N, TZHBCIMMSKGBFX-QJUCYEIZSA-N | commercial name, resolvers disagree | name-only |
| `4-pregnene-11β,17,21-triol-3,20-dione 21-acetate` | FAILED | - | commercial name names two things |
| `4-pregnene-11β,17,21-triol-3,20-dione 21-hemisuccinate` | AFLWPAGYTPJSEY-CODXZCKSSA-N | commercial name only | name-only |
| `4-pregnene-17-ol-3,20-dione caproate` | DOMWKUIIPQCAJU-LJHIYBGHSA-N | commercial name only | name-only |
| `4-pregnene-17α,21-diol-3,11,20-trione` | MFYSYFVPBJMHGN-ZPOLXVRWSA-N | commercial name only | name-only |
| `4-pregnene-17α-ol-3,20-dione` | DBPWSSGDRRHUNT-CEGNMAFCSA-N | commercial name only | name-only |
| `4-pregnene-3,20-dione` | RJKFOVLPORLFTN-LEKSSAKUSA-N | PubChem+CIR |  |
| `4-pregnene-6α-methyl-17-ol-3,20-dione` | FRQMUZJSZHZSGN-HBNHAYAOSA-N | commercial name only | name-only |
| `4-pregnene-6α-methyl-17-ol-3,20-dione acetate` | PSGAAPLEWMOORI-PEINSRQWSA-N | commercial name only | name-only |
| `5-androstene-3β-ol-17-one glucuronide` | GLONBVCUAVPJFV-UHFFFAOYSA-N | commercial name only | name-only |
| `5-androstene-3β-ol-17-one sulfate` | FAILED | - | commercial name names two things |
| `5α-androstane-17α-methyl-17β-ol-3,2c-pyrazole` | LKAJKIOFIWVMDJ-IYRCEVNGSA-N | commercial name only | name-only |
| `5α-androstane-17α-methyl-17β-ol-3one` | WYZDXEKUWRCKOB-YDSAWKJFSA-N | commercial name only | name-only |
| `5α-androstane-17β-ol-3-one` | FAILED | - | commercial name names two things |
| `5α-androstane-1α-methyl-17β-ol-3one` | UXYRZJKIQKRJCF-TZPFWLJSSA-N | commercial name only | name-only |
| `5α-androstane-3,17-dione` | RAJWOBJTTGJROA-WZNAKSSCSA-N | PubChem+CIR |  |
| `5α-androstane-3α-ol-17-one sulfate` | ZMITXKRGXGRMKS-HLUDHZFRSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity ZMITXKRGXGRMKS |
| `5α-androstane-3β-ol-17-one glucuronide` | VFUIRAVTUVCQTF-PALHZPRPSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity VFUIRAVTUVCQTF |
| `5α-androstane-3β-ol-17-one sulfate` | ZMITXKRGXGRMKS-LUJOEAJASA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity ZMITXKRGXGRMKS |
| `5α-estrane-3α-ol-17-one glucuronide` | JGNAYBFRYHOLMC-UHFFFAOYSA-N | commercial name only | name-only; stereo-pair: 2 names share connectivity JGNAYBFRYHOLMC |
| `5β-androstane-3α-ol-17-one glucuronide` | VFUIRAVTUVCQTF-SDHZCXLISA-N, VFUIRAVTUVCQTF-UHFFFAOYSA-N | commercial name, resolvers disagree | name-only; stereo; stereo-pair: 2 names share connectivity VFUIRAVTUVCQTF |
| `5β-estrane-3α-ol-17-one glucuronide` | JGNAYBFRYHOLMC-HAVPAZTRSA-N, JGNAYBFRYHOLMC-UHFFFAOYSA-N | commercial name, resolvers disagree | name-only; stereo; stereo-pair: 2 names share connectivity JGNAYBFRYHOLMC |
| `pregna-1,4-diene-3,20-dione, 9-fluoro-11,21-dihydroxy-16,17-((1-methylethylidene)bis(oxy))-, (11β,16α)-` | YNDXUCZADRHECN-JNQJZLCISA-N | commercial name only | name-only |

## What "resolved clean" means, and what NEEDS CONFIRMATION means

A name is RESOLVED only where a resolver returned exactly one InChIKey FOR THE SYSTEMATIC
NAME. Resolving from the commercial name instead is not a resolution and is reported as
ambiguous, because a commercial name is the weaker identifier and several here name two
compounds at once.

**NEEDS CONFIRMATION** marks a clean resolution that shares connectivity with another name
in this corpus, differing only in stereochemistry. In these the resolvers DID distinguish
the pair, returning different stereo blocks, so the answer looks right - and that is the
point. It is the one case where a WRONG answer would also look right, because the only thing
separating the two compounds is a stereo block nothing downstream inspects, and the two have
different cross sections. These do not need re-resolving. They need a chemist to confirm
which stereoisomer is which, once, and record it.

REPRODUCIBILITY. This report is a SNAPSHOT of two live public services. Unlike the model's
own digests, re-running it on another day may give different answers if PubChem or CIR change
what they return. That is a further reason no assignment here belongs in the corpus without a
person putting it there.

## Flags

- **name-only** - 67
- **stereo-pair** - 18
- **stereo** - 10
- **commercial name names two things** - 8
