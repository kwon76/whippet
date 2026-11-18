# Claude Code 작업 지시서 v2: Whippet 리비전 (first fold 구현 → end-to-end 실험 → 논문 반영)

이 문서 전체를 Claude Code에 붙여 넣으면 된다. 작업 디렉터리는 이 저장소의 루트다.
`CLAUDE_CODE_PROMPT.md`(v1)의 내용을 대체한다. v1과 다른 점은 다음과 같다.
- 리뷰어 질문 Q1–Q5를 기준으로 실험을 설계한다.
- 작업을 단계별로 나누고 각 단계에 통과 조건(gate)을 둔다.
- 코드를 직접 읽고 확인한 구현 세부를 반영했다.
- 결과가 논문 주장을 흔들 때의 결정 절차를 정했다.
- 논문을 고칠 위치를 문장 단위로 지정했다.

---

## 0. 배경

이 저장소는 CRYPTO 투고용 논문 "Whippet: A Shorter Tail for Lattice-Based Polynomial Commitments"의 리비전 작업 공간이다.
- `tex/`: 논문 소스
- `model/`: 증명 크기 모델
- `proto/`: Python 참조 구현. tail만 구현되어 있다.
- `build.py`: JSON에서 수치·표를 만들고 `out/whippet.pdf`를 빌드한다.

1차 리뷰 판정은 borderline / weak reject였다. 남은 핵심 지적은 두 가지다(`NOTES.md`의 "Revision after the first review" 참고).
1. 구현은 tail만 다루고, 입력은 first fold 출력과 길이·노름 한계만 같은 합성 witness다. 실제 다항식에서 시작한 commit → open → verify는 검증된 적이 없다.
2. 실제 witness에서 정직한 prover의 시도 횟수, 실패율, 시간·메모리 자료가 없다.

이번 리비전의 목표는 프로토콜이나 파라미터를 바꾸는 것이 아니다. **first fold를 구현하고, 실제 다항식으로 PCS 전체를 실행하고, 그 결과로 논문의 주장을 확인하거나 고치는 것**이다.

## 1. 리뷰어 질문과 논문 반영 위치

모든 실험은 아래 질문 중 적어도 하나에 답해야 한다. 질문에 답하지 않는 실험은 하지 마라.

| | 질문 | 답이 들어갈 곳 (현재 문장) |
|---|---|---|
| Q1 | 실제 다항식에서 시작한 PCS 전체가 정확하고 안전하게 동작하는가? (정직한 증명은 통과, 변조는 모두 거부) | 초록 "…the first fold is not implemented.", 1장 *Validation* 문단, 1.2 *Limitations* "The first fold is modelled but not implemented…", 7장 *Scope* |
| Q2 | 구현이 쓴 바이트(first fold 몫 포함)가 모델과 얼마나 다르고, 왜 다른가? | 7장 *Sizes*, `tab:impl` 캡션 "The first fold would add the commitment v and 12 bytes." |
| Q3 | 최적화(증명 기반/483) 설정에서 실제 first fold witness의 노름은 모델의 기대 노름 한계 대비 얼마인가? 어떤 입력에서 실패하며, 빈도와 영향은? | 6.6 "The optimised parameters" 문단, 7장 *Attempts and norms*, 8장 "Completeness and real witnesses" |
| Q4 | 보장 설정에서 최악의 다항식을 넣어도 정리 6.3대로 동작하는가? (모든 레벨에서 `fill ≤ 1`, 기대 시도 ≤ 2) 실제 증명은 상한보다 얼마나 짧은가? | 정리 5.2 뒤 "…whatever the polynomial.", 6.6 *Cost* 마지막 문장, 7장 |
| Q5 | 참조 구현 기준으로 commit / open / verify의 시간과 메모리는 얼마인가? | 7장(새 표), 1.2 "we have not measured verification times" |

마지막 보고에는 Q1–Q5 각각에 대해 수치가 들어간 한두 문장의 답이 있어야 한다.

## 2. 규칙

v1의 규칙을 모두 유지한다.
- 수치를 지어내지 않는다.
- 논문과 코드가 다르면 `NOTES.md`에 기록하고 근거를 들어 보고한다.
- 검증기를 약하게 만들지 않는다.
- 덮어쓰기 전에 `results_prev/`에 백업한다.
- 영국식 철자를 쓰고, 과장 표현을 쓰지 않으며, 수치는 매크로로만 넣는다.
- 모호한 명세는 선택한 해석과 근거를 `NOTES.md`에 적는다.
- 정리 4.2, 5.2, 6.3과 어긋나는 선택은 하지 않는다.

추가 규칙:
1. **예측을 먼저 기록한다.** 각 실험 묶음을 돌리기 전에 `NOTES.md`에 "예측: …"을 한 줄씩 적는다. 결과가 나오면 예측과 비교한다. 예를 들어 최적화 설정에서 모든 블록이 같은 다항식은 `z = (Σ c_i) s_1`이 되어 tail 레벨 1의 투영이 반복 실패할 것으로 예상할 수 있다. 사후 해석을 막기 위한 규칙이다.
2. **통계를 정직하게 보고한다.** 실패가 0/20이면 "실패 없음"이라고 쓰지 않는다. "0/20 (95% 상한 14%)"처럼 Clopper–Pearson 상한을 함께 쓴다. 평균 시도 횟수에는 최댓값과 실행 수를 같이 쓴다.
3. **결정 게이트를 지킨다.** 결과가 논문의 헤드라인 주장을 흔들면 초록·1장의 헤드라인 문장을 고치기 전에 멈추고 보고한다. 헤드라인은 "최적화 설정 40.7 KB가 실제 다항식에서도 동작한다"는 전제다. 예를 들어 어떤 입력 종류에서든 실패율이 5%를 넘는 경우가 해당한다. 보고에는 선택지 2–3개와 각각의 크기 비용을 넣는다(§6.4 참고). 7장과 8장의 사실 기술은 게이트 없이 고쳐도 된다.
4. **기존 동작을 보존한다.** `tail.py`를 일반화할 때 기존 호출 경로의 출력이 비트 단위로 같아야 한다. `run_tests.py`를 다시 돌려 `tests.json`의 크기와 검증 결과가 바이트까지 같은지 확인한다. 다르면 이유를 `NOTES.md`에 쓴다.
5. **재현 정보를 남긴다.** `proto/e2e.json`에 다음을 기록한다: git 커밋 해시, CPU, 코어 수, RAM, Python·NumPy 버전, 시드, 각 실행의 시작·종료 시각. 큰 실험은 백그라운드로 돌리고, 실행 단위로 JSON에 저장해 중단돼도 이어갈 수 있게 한다.
6. **컨텍스트 압축에 대비한다.** 각 단계가 끝날 때마다 `NOTES.md`의 "State on <날짜>" 절에 현재 상태, 다음 할 일, 열린 문제를 적는다.
7. 커밋은 단계가 끝날 때마다 작게 나눠서 한다. 브랜치 규칙은 세션 지시를 따른다.

## 3. 단계 0 — 준비 (게이트: 현재 상태로 빌드와 테스트 통과)

1. `tex/`가 저장소 루트에 있는지 확인한다. `build.py`는 `ROOT/tex`를 읽는다. 없으면 사용자가 올린 `tex.zip`을 루트에 푼다(`tex/00_preamble.tex` … `tex/91_appendix.tex`).
2. `pdflatex`를 확인한다. 없으면 TikZ, pgfplots, cleveref, enumitem, booktabs가 들어 있는 TeX Live 패키지를 설치한다. 설치할 수 없으면 그렇다고 보고하고 빌드 확인은 마지막에 한다.
3. 다음 파일을 읽는다.
   - 문서: `README.md`, `NOTES.md`
   - 논문: `tex/04_tail.tex`, `tex/05_pcs.tex`(Figure `fig:greyhound`, 식 `eq:eval`, `eq:gh-system`, `eq:gh-open`, 정리 `thm:pcs`), `tex/06_params.tex`(6.6 `sec:completeness`), `tex/07_impl.tex`, `tex/08_discussion.tex`, `tex/01_intro.tex`
   - 코드: `proto/tail.py`, `proto/ring.py`, `model/gh.py`, `model/wmodel.py`(`WOpt`, `search`, `_wlevel_gc`), `proto/run_full.py`, `proto/run_retries.py`
4. 기존 결과 파일을 `results_prev/`에 복사한다: `model/final.json`, `proto/full_*.json`, `proto/retries.json`, `proto/tests.json`, `proto/pairs.json`.
5. `cd proto && python3 run_tests.py`와 `python3 build.py`(`--draft` 없이)를 현재 상태 그대로 돌린다. 이 PDF가 비교 기준(baseline)이다. Overfull 경고와 참조 경고를 기록해 둔다.
6. 기계 사양을 기록한다. 직전 세션의 컨테이너는 4코어, 15 GB RAM이었다.

## 4. 단계 1 — first fold 구현 (`proto/pcs.py`)

### 4.1 명세

5장의 Figure `fig:greyhound`와 식 (`eq:eval`), (`eq:gh-system`)가 명세다.

- **기호**: `d = 64`, `q = 2^32 - 99`, `N = m*r*d`. 블록 `f_i`는 환 원소 `m`개다.
- **계수 표현**: `[0, q)` 계수를 중심화 대표 `(-q/2, q/2]`로 바꾼 뒤 자릿수로 분해한다. 모델의 최상위 자릿수 한계 `T_0 = (q/2)/b_0^(δ0-1) + 1`이 이 규약을 전제한다. `NOTES.md`에 명시한다.
- **자릿수 분해**: `s_i = G^{-1}_{b0,m}(f_i)`. 마지막 자릿수는 남은 몫이다(`tail.digits`와 같은 규약).
- **커밋과 평가**: `t_i, t^H_i, l_i, t̂_i, u, e_u, w_i, ŵ, v, e_v, ŷ, z, e`는 v1의 정의를 따른다. 평가 쪽 정의는 다음과 같다.
  - `a = (x^{dj})_{j<m}`
  - `b = (x^{dmi})_{i<r}`
  - `ψ_x = σ_{-1}(Σ_k x^k X^k)`
  - `ā = G^T a`

### 4.2 `LevelStatement` 재사용 — 코드에서 확인한 대응

`tail.LevelStatement.agg`를 읽고 확인한 내용이다. 구현 전에 직접 다시 확인하라.

- **witness 순서**: `[z0 | z1 | t̂ | ŵ | e0 | e1 | d1 | d2]`이다. `Bm`(t̂의 외부 키)이 `B`이고 `Cm`(ŵ의 외부 키)이 `D`이므로 `d1 = e_u`, `d2 = e_v`다. 논문 Figure의 나열 순서 `(ŵ, t̂, z0, z1, e_v, e_u, e0, e1)`와는 다르다. 수학적으로는 무관하지만 코드는 `agg`의 순서를 따라야 한다.
- **`lv` 대응**:

  | 키 | 값 |
  |---|---|
  | `n` | `δ0*m` |
  | `r` | `r` |
  | `kap` | 내부 랭크 `n` |
  | `k1` | `n1` |
  | `t1` | `kt` |
  | `tw` | `kw` |
  | `b`, `b1`, `bw` | `2^lb` |
  | `di` | `Dt` (= ρ_t) |
  | `do` | `Du` (= ρ_u) |
  | `dec` | `True` |
  | `wmode` | `'commit'` |
  | `m2` | witness 길이 |

  값은 `gh.instance` 또는 `gh.instance_gc`의 반환값에서 가져온다.
- **`E`**: `ā`를 `F` 원소로 넣는다(좌표 0만 채우고 나머지는 0).
- **`PhiP`**: `PhiP_i = ψ_x · b_i ∈ R_q`다. tail에서도 `PhiP`는 이미 환 원소 배열 `(r, 64)`이다(`tail.py`의 `PhiP = Pt[:, 0]`). 따라서 타입을 일반화할 필요는 없을 것으로 보인다. 확인하라.
- **`claim`**: `ŷ ∈ F` 전체다. `τ_0(ŷ) = y`는 verifier가 채운다.
- **ŵ 자릿수 규칙**: 최적화 설정의 `gh.first_fold`은 table rule(`widths(LQ, True)`: 마지막 자릿수가 더 넓을 수 있음)을 쓴다. 보장 설정의 `gh.first_fold_gc`는 `table_rule` 인자를 무시하고 ceil 규칙을 쓴다. `tail.digits(x, base, k)`는 마지막 자릿수가 남은 몫이므로 `k = kw`만 맞으면 두 규칙을 모두 처리할 것이다. 노름 한계가 모델과 일치하는지 확인하라.
- **제약식**: 조밀 벡터로 하나씩 만들지 말고 `agg`가 만드는 집계 벡터 하나만 만든다.

### 4.3 키, Fiat–Shamir, 인터페이스

- **키**: `A, B, D`는 `tail.crs`로 시드에서 유도하고, 태그는 tail 레벨과 겹치지 않게 정한다(예: `b"pcs|A"`).
- **Fiat–Shamir 흡수 순서**:
  1. 공개 파라미터(모양, 랭크, 밑, `Dt`, `Du`, 모든 tail 레벨의 schedule. `bind_statement`의 파라미터 부분과 같은 형식)
  2. `u`, `x`, `y`
  3. 그다음 `v`, `ŷ`의 3좌표, first fold nonce
  4. `c_1..c_r` 유도
  5. 같은 transcript로 tail 진입

  `tail.prove`와 `tail.verify`에 "이미 묶인 transcript를 받는" 선택 인자를 추가한다. 기본값일 때의 동작은 기존과 비트 단위로 같아야 한다.
- **설정별 first fold 동작**:
  - 최적화 설정: 모델이 first fold에서 챌린지 재시도를 가정하는지 `gh.first_fold`과 `wmodel`에서 확인하고, 그대로 구현한다. 재시도가 없다고 가정한다면 first fold nonce는 0으로 고정되고, 노름 초과는 tail 레벨 1의 투영 재시도로 나타난다. 이 해석을 `NOTES.md`에 적는다.
  - 보장 설정: `‖z‖² + ‖e‖² ≤ g['qmax']`가 될 때까지 챌린지를 다시 뽑는다. 6.6절 마지막 문단의 내용이다.
  - nonce 상한은 tail과 같은 250으로 둔다. 초과하면 예외를 던지고 레벨을 기록한다. 상한을 바꾸지 마라.
- **인터페이스**: `setup(params)`, `commit(f) -> (u, state)`, `open(state, x) -> (y, proof_bytes)`, `verify(u, x, y, proof_bytes) -> bool`.
- **증명 바이트**: 증명에는 `v`, `ŷ`의 3좌표, nonce, tail 증명이 들어간다. `Writer.sizes`에 `first_fold` 항목을 추가한다.
- **바이트 회계**: `gh.own_bits`는 챌린지 시드 2개(256비트)를 first fold 몫에 포함한다. 반면 논문의 `tab:impl`은 "without challenge seeds"다. 구현 바이트와 모델을 비교할 때 무엇을 포함하는지 한 가지 규약으로 통일하고 표 캡션에 쓴다.

### 4.4 작은 크기용 모양

`gh.SHAPE`의 튜플은 `(m, r, lb0, δ0, lb, dl)`이다. 작은 크기 `N = 2^18, 2^20, 2^22, 2^24`를 추가한다.
- `m*r ≥ N/64`, `m/r ≈ 9.5`
- `lb0 = 4`(`b0 = 16`, `δ0 = 8`), `lb = 6`
- 정확히 나누어떨어지지 않으면 0으로 채운다(padding). 선택한 `(m, r)`과 padding 양을 기록한다.

각 크기의 tail schedule은 `W.set_bkz(483)` 뒤에 아래처럼 `W.search`를 직접 불러 구하고 `model/small_sched.json`에 캐시한다. `final.json`에는 쓰지 않는다.
- 최적화: `gh.instance(logN, True, 'A')` + `W.WOpt(acc='A')`
- 보장: `gh.instance_gc(logN, 2.0, 128**-0.5, 4, 6)` + `W.WOpt(acc='A', gc=True)`
- 중간 옵션 P: `W.WOpt(acc='A', Pm=True)`. 6.6절 마지막 문단과 `final.json`의 `whip_gh30_A_483_P`에 해당한다. 단계 2에서 최적화 설정이 실패할 때만 필요하다.

랭크 고정점이 수렴하지 않거나 `instance_gc`가 `None`을 반환하면, 그 원인(어느 노름이 `q`를 넘는지)을 보고한다.

### 4.5 게이트 1 (`N = 2^18`, 최적화·보장 두 설정)

- 정직한 증명이 통과한다. 균등 입력과 모든 블록이 같은 입력 각각 3회 이상.
- §5.3의 음성 테스트가 모두 거부된다.
- 구현 바이트와 모델(`gh.own_bits` + tail 모델)의 차이가 기존 tail 수준(0.1% 내외)이고, 차이의 원인이 항목별로 설명된다.
- `run_tests.py`가 기존과 동일하게 통과한다(규칙 4).
- 커밋한다.

## 5. 단계 2 — 실험 (`proto/run_e2e.py`, 결과는 `proto/e2e.json`)

### 5.1 조합

**크기와 실행 횟수**
- `N = 2^18, 2^20, 2^22`: 각 조합 20회
- `N = 2^24`: 각 조합 5회 이상
- `N = 2^26`: 가능하면. 단계 3의 성능 작업 결과에 달려 있다.
- 실행마다 다항식 시드와 평가점을 바꾼다.

**설정**
- 최적화/483
- 보장/483
- 최적화가 어떤 입력에서 실패하면 중간 옵션 P/483을 추가한다.

**입력 종류**

| 기호 | 입력 | 확인하려는 것 |
|---|---|---|
| U | `[0,q)` 균등 | 모델의 기대 노름 가정에 유리한 경우 |
| S | 희소(1%만 0이 아님) | 자릿수 대부분이 0 |
| E | 모든 계수가 같은 값 | 모든 블록이 같아 `z = (Σc_i)s_1` |
| P7 | 주기 7로 반복 | 블록 간 상관이 약함 |
| Pd | 주기 `d·m`로 반복 | 모든 블록이 같음(E의 일반형) |
| B8 | 8비트 작은 계수 | 상위 자릿수가 0 |
| Mx | 모든 자릿수가 최댓값(부호 무작위) | 보장 한계 `B_s`에 닿음 |
| MxE | Mx와 같은 크기이고 모든 블록이 같음 | 보장 설정에서 `‖z‖`의 최악 |

**평가점**: 무작위 `x`에 더해 `x = 0, 1, q-1`을 각 1회 넣는다. 구조가 단순한 statement에서도 깨지지 않는지 확인하기 위해서다.

### 5.2 기록 항목

| 항목 | 내용 |
|---|---|
| 바이트 | 전체, `first_fold`, 레벨별. 모델 값과 차이(%) |
| 노름 | first fold 출력(레벨 1 입력)의 `‖w‖/β_1`(모델 `g['beta']` 대비)과 `‖(z,e)‖/γ̄`. 레벨별 `‖w_j‖/β_j`(`tail.STATS`의 `fill`) |
| 시도 | first fold 챌린지, 레벨별 `jl_tries`, `c_tries`. 평균·최대·히스토그램 |
| 실패 | 예외 종류와 레벨. 실패율은 Clopper–Pearson 95% 구간과 함께 |
| 자원 | commit / open / verify 각각의 wall time과 최대 RSS(`resource.getrusage`). 단계마다 별도 프로세스로 재거나, 단계 사이의 RSS 증가분을 기록 |

### 5.3 음성 테스트 (각 크기·설정에서 최소 1회, 모두 거부되어야 함)

- 틀린 `y`(`y+1`, 그리고 비정규 표현 `y+q`가 입력으로 가능하면 그것도)
- 틀린 `x`
- 다른 다항식의 `u`
- 같은 증명을 다른 `(x, y)`에 재사용
- 공개 파라미터 한 개 변경
- `ŷ` 좌표 변조
- 임의 비트 반전 60회, 절단, 1바이트 덧붙이기

### 5.4 기존 tail 실험 재실행 (최종 코드 기준)

- `run_tests.py`
- `run_retries.py retries.json 20`
- `README.md`의 `run_full.py` 다섯 줄. 각 실행은 약 10분, 3 GB.
- 가능하면 `whip_gh30_A_483_full`을 20회 반복해 시도 횟수 분포를 구하고 `proto/full_Q_rep.json`에 저장한다.

### 5.5 게이트 2

`NOTES.md`에 다음을 정리한다.
- 결과 요약 표
- Q1–Q4에 대한 수치가 있는 답
- 예측과 결과의 대조

최적화 설정에서 실패가 있으면 여기서 규칙 3의 결정 게이트를 적용한다.

## 6. 단계 3 — 성능 (큰 크기를 위해. 필요할 때만)

1. **먼저 측정한다.** `N = 2^18, 2^20`의 commit / open / verify 시간과 RSS로 `2^22–2^26`을 외삽해 `NOTES.md`에 적는다. `2^24`의 commit이 1시간을 넘을 것으로 보이면 이 단계를 단계 2의 `2^24` 실행보다 먼저 한다.
2. **곱셈을 가속한다.** 대상은 `ring.matvec_small`, `mul_small`, `mul`이다. 선택지는 다음과 같다.
   - numba 또는 C 확장(ctypes/cffi)
   - 정확한 정수 행렬곱 기법: 균등 행렬을 16비트 limb로 쪼개고 negacyclic 곱을 Toeplitz 전개해 int64/float64 BLAS로 처리한다. 누적합이 `2^53`(float64) 또는 `2^63`(int64)을 넘지 않음을 식으로 보인다.
   - 보조 NTT 소수와 CRT. `q = 2^32-99`에서는 `X^64+1`이 완전 분해되지 않기 때문이다.

   어느 방법이든 무작위 입력 1000개 이상에서 기존 `ring.*`와 비트 단위로 같아야 한다. 이 테스트를 `run_tests.py`에 추가한다.
3. **메모리를 줄인다.** `s_i`를 블록 단위로 스트리밍한다. 커밋은 블록별로 하고, `z`는 `z += c_i s_i`로 누적하며, 자릿수는 int8로 두거나 필요할 때 재계산한다. 15 GB 기계에서 `2^28`은 시도할 만하다. `2^30`은 다항식만 4 GB이므로 스트리밍이 끝난 뒤에만 시도한다.
4. **실제로 돌린 크기만 보고한다.** 같은 기계에서 잰 값이 없으면 Greyhound·Akita와 시간을 비교하는 문장을 쓰지 않는다.

## 7. 단계 4 — 논문 반영

### 7.1 생성 경로

- `build.py`가 `proto/e2e.json`을 읽어 매크로(`\nEe...`)와 표(`gen/tab_e2e_*.tex`)를 만든다.
- 최종 빌드는 `--draft` 없이 한다. `--draft`는 없는 결과를 다른 결과로 채우는 자리표시자 모드다.
- `e2e.json`이 없거나 필요한 키가 없으면 조용히 넘어가지 말고 빌드를 실패시킨다.

### 7.2 고칠 곳

| 위치 | 현재 | 고칠 내용 |
|---|---|---|
| 초록 | "The model agrees to within 0.1% with a reference implementation of the recursion on synthetic witnesses of full size; the first fold is not implemented." | end-to-end 실행 범위(크기, 설정, 입력 수)와 모델 대비 편차. 돌리지 못한 크기를 명시 |
| 1장 *Validation* | 합성 witness만 언급 | 실제 다항식 실행 결과를 한두 문장으로 |
| 1.2 *Limitations* | "The first fold is modelled but not implemented … main gap …" | 남은 간극을 사실대로(예: `2^28`/`2^30` end-to-end 미실행, 최적화 설정의 구조적 입력 실패율, Python 시간의 의미). "we have not measured verification times"는 "참조 구현의 시간만 쟀고 비교는 하지 않는다"는 취지로 |
| 5.2 정리 `thm:pcs` 뒤 | "…whatever the polynomial." | 이론 문장은 그대로 두고, 최악 입력(Mx, MxE) 실험을 7장 참조로 덧붙인다 |
| 6.6 *Cost* 마지막 문장 | "proofs for witnesses below the bound are shorter (Section 7)" | 실제 다항식에서 상한 대비 몇 % 짧은지 |
| 6.6 *The optimised parameters* | 무작위·구조적 합성 witness만 언급 | 실제 first fold witness 결과(Q3). P 옵션을 실험했다면 그 결과도 |
| 7장 *Scope* | "The first fold is not implemented…" | 실제로 한 일로 다시 쓴다 |
| 7장 새 문단 *End-to-end runs* + 표 | 없음 | (a) `N`별 바이트, 모델, 편차, first fold 몫 (b) 입력 종류별 성공 수/실행 수, 평균·최대 시도 (c) 레벨별 `‖w‖/β`와 시도 (d) 시간·RSS |
| `tab:impl` 캡션 | "The first fold would add the commitment v and 12 bytes." | 실측 first fold 몫으로 교체하거나 새 표를 참조 |
| 8장 "Completeness and real witnesses" | "…has not been tested; an implementation of the first fold is needed…" | 결과로 다시 쓴다. 남은 열린 문제(더 촘촘한 분산 기반 한계 등)는 남긴다 |
| `README.md` "Status of the numbers" | "covers the tail (not the first fold)" | 갱신 |

### 7.3 지면

CRYPTO 분량을 고려해 본문 7장에는 표 두 개까지만 둔다(크기 대 모델, 입력 종류별 결과). 레벨별 노름·시도와 시간·메모리 표는 부록으로 보내고 본문에서 참조한다.

### 7.4 결과가 주장과 어긋날 때 (규칙 3)

최적화 설정이 구조적 다항식(E, Pd, MxE)에서 자주 실패하면, 아래 선택지를 크기 비용과 함께 보고하고 사용자 결정을 기다린다. 결정 없이 헤드라인을 바꾸지 마라.
1. 헤드라인은 유지하되, 최적화 설정의 완전성이 휴리스틱이며 어떤 입력에서 실패율이 얼마인지 초록과 1.2에 명시한다.
2. 보장 설정(47.4 KB)을 "안전한 기본값"으로 앞세운다.
3. 중간 옵션 P(`\nAblP` KB 추가)를 실제 다항식에서 검증된 절충안으로 제시한다.

### 7.5 빌드 확인

`python3 build.py`로 다시 빌드한다. 다음을 확인한다.
- `Overfull`이 4pt 이하다.
- 미정의 참조·인용 경고가 없다.
- 단계 0의 baseline PDF와 비교해 의도하지 않은 수치 변화가 없다. 기존 매크로 값이 바뀌었다면 이유를 적는다.

### 7.6 리뷰 대응 자료

- `out/changes.md`: 고친 문장마다 "이전 → 이후 → 근거(어느 JSON의 어느 키)"를 적는다.
- `out/response_points.md`: 리뷰의 각 지적에 대해 무엇을 했는지와 논문 위치를 적는다. 저장소에 리뷰 원문이 있으면 그것을 쓰고, 없으면 `NOTES.md`의 요약을 쓴다. 대응 편지 문체가 아니라 사실 목록으로 쓴다.

## 8. 완료 기준

1. `proto/pcs.py`의 commit / open / verify가 실제로 돌린 모든 크기에서 통과하고, 음성 테스트가 모두 거부된다.
2. `proto/e2e.json`과 요약 표가 있다. 모델 대비 바이트 차이와 그 원인이 항목별로 설명되어 있다.
3. 기존 테스트와 재실행 결과가 최종 코드 기준으로 갱신되어 있고, 이전 결과는 `results_prev/`에 있다.
4. 논문이 새 결과로 빌드되고, 돌리지 못한 항목이 *Limitations*에 남아 있다.
5. `NOTES.md`에 다음이 정리되어 있다: 한 일, 논문과 코드의 차이, 예측과 결과, 남은 일.

## 9. 최종 보고 (한국어)

- 실제로 돌린 `(N, 설정, 입력)` 목록과 각 실행 수
- Q1–Q5 각각에 대한 한두 문장의 답(수치 포함)
- 실패 횟수와 레벨(Clopper–Pearson 상한 포함)
- 모델 대비 최대 편차와 원인
- 최적화 설정에서 관찰된 `‖w‖/β`의 범위, 최대 시도 횟수, 실제 first fold witness 노름 대비 모델 기대 노름 한계
- 결정 게이트에 걸린 항목과 제안한 선택지
- 논문에서 고친 문장 목록(`out/changes.md` 요약)
- 돌리지 못한 것과 그 이유

## 10. 알려진 주의점

- `tail.random_instance`는 반올림 때문에 한계를 1–2% 넘는 witness를 만들 수 있다. 합성 witness가 필요하면 `run_full.py`나 `run_retries.py`의 생성기를 쓴다. e2e 실험은 실제 다항식이므로 이 문제가 없다.
- 보장 설정에서 `fill > 1`이 한 번이라도 나오면 구현이나 한계식(6.6절 식 (13)–(14), `wmodel._wlevel_gc`, `gh.first_fold_gc`) 중 하나가 틀린 것이다. 원인을 찾을 때까지 다음 단계로 가지 않는다.
- `model/final.json`의 `trace`가 schedule이다. `bjl`은 투영 한계, `qmax`는 `(z, e)`의 결합 한계, `joint`는 종단 레벨에서 결합 검사를 쓴다는 표시다.
- `model/run_final.py` 전체 재계산은 2코어에서 수 시간 걸린다. 필요한 항목만 `W.search`로 따로 계산한다.
- `NOTES.md` 끝의 Akita 수치 "CHECK"(66.1 vs 67.2)는 이번 작업 범위가 아니다. 원문을 확인할 수 없으면 건드리지 않는다.
