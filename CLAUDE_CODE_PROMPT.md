# Claude Code 작업 지시서: Whippet 논문의 end-to-end 구현과 실험

아래 전체를 Claude Code에 그대로 붙여 넣으면 됩니다. 작업 디렉터리는 `whippet_code.zip`을 푼 폴더입니다.

---

## 0. 배경과 목표

이 저장소는 CRYPTO 투고를 준비 중인 논문 "Whippet: A Shorter Tail for Lattice-Based Polynomial Commitments"의 소스(`tex/`), 증명 크기 모델(`model/`), Python 참조 구현(`proto/`), 빌드 스크립트(`build.py`)다. 먼저 `README.md`, `NOTES.md`, `tex/04_tail.tex`, `tex/05_pcs.tex`, `tex/06_params.tex`의 6.6절, `tex/07_impl.tex`, `proto/tail.py`, `model/gh.py`, `model/wmodel.py`를 읽어라.

리뷰어의 핵심 지적은 두 가지다.

1. 참조 구현은 재귀 부분(tail)만 구현하고, 입력은 "first fold가 만드는 것과 길이·노름 한계가 같은 합성 witness"다. **실제 다항식에서 출발한 polynomial commitment 전체 실행은 검증된 적이 없다.**
2. 정직한 prover가 실제 witness에서 몇 번의 재시도로 끝나는지, 실패(nonce 소진)가 있는지, 시간·메모리는 얼마인지에 대한 자료가 없다.

이번 작업의 목표는 **first fold를 구현해 tail과 연결하고, 실제 다항식으로 commit → open → verify 전체를 돌려 그 결과를 논문 7절에 반영하는 것**이다. 프로토콜이나 파라미터를 바꾸는 작업이 아니다.

## 1. 절대 규칙

- 숫자를 지어내지 마라. 논문에 들어가는 모든 수치는 JSON 결과 파일에서 `build.py`가 생성해야 한다. 돌리지 못한 크기는 "돌리지 못했다"고 쓴다.
- 논문과 코드가 다르면 조용히 한쪽을 고치지 말고, 차이를 `NOTES.md`에 기록한 뒤 어느 쪽이 맞는지 근거와 함께 보고하고 양쪽을 일치시켜라.
- 검증기를 약하게 만들지 마라. 기존 테스트(`proto/run_tests.py`)는 계속 통과해야 한다.
- 기존 결과 파일(`model/final.json`, `proto/full_*.json`, `proto/retries.json`, `proto/tests.json`)을 덮어쓰기 전에 `results_prev/`에 복사해 둬라.
- 논문 문체: 영국식 철자, 과장 표현 금지, 수치는 매크로(`\n...`)로만. 한계는 숨기지 말고 본문에 쓴다.

- 명세가 모호한 곳은 추측으로 넘어가지 말고, 선택한 해석과 근거를 `NOTES.md`에 적은 뒤 진행하라. 논문의 정리(4.2, 5.2, 6.3)와 어긋나는 선택은 하지 마라.
- 실행 환경: Python 3.10 이상, NumPy. 큰 실험은 백그라운드로 돌리고 중간 결과를 JSON에 계속 저장해 중단돼도 이어갈 수 있게 하라.

## 2. 작업 A — first fold 구현 (`proto/pcs.py` 신규)

논문 5절(Figure 4, 식 (10)–(12))이 명세다. 요지:

- 다항식 `f`는 `Z_q` 계수 `N = m*r*d`개(`d = 64`, `q = 2^32 - 99`). 블록 `f_i`(환 원소 `m`개)를 밑 `b_0`의 중심화 자릿수 `delta_0`개로 분해해 `s_i`(환 원소 `delta_0*m`개)를 만든다. 마지막 자릿수는 "남은 몫"이다(`tail.py`의 `digits`와 같은 규약).
- 커밋: `t_i = A s_i`, `t_i^H = hi_{rho_t}(t_i)`, `l_i = lo_{rho_t}(t_i)`, `t_hat_i = digits(t_i^H)`, `u = hi_{rho_u}(B t_hat)`, `e_u = lo_{rho_u}(B t_hat)`. 커밋먼트는 `u`.
- 평가 `f(x) = y`: `a = (x^{d j})_{j<m}`, `b = (x^{d m i})_{i<r}`, `psi_x = sigma_{-1}(sum_k x^k X^k)`, `a_bar = G^T a`. `w_i = <a_bar, s_i>`, `w_hat = digits(w)`, `v = hi_{rho_u}(D w_hat)`, `e_v = lo(...)`, `y_hat = Theta(psi_x * sum_i b_i w_i)`(좌표 0은 `y`이므로 나머지 3개만 전송). 챌린지 `c_1..c_r`, `z = sum c_i s_i`, `e = sum c_i l_i`, 밑 `b`로 두 부분 분해.
- tail에 넘기는 관계는 식 (11)이고 이것은 식 (7)과 이름만 다르다(`w_hat <-> h_hat`, `e_u <-> d_1`, `e_v <-> d_2`, `a_bar <-> E`, `psi_x*b_i <-> Phi'_i`, `y_hat <-> v`). 따라서 **`tail.py`의 `LevelStatement`를 재사용할 수 있는지 먼저 확인하라.** witness 조각의 순서는 `LevelStatement.agg`가 기대하는 순서 `[z0 | z1 | t_hat | w_hat | e 부분 | d_1 | d_2]`와 정확히 맞춰야 한다. `a_bar`의 성분은 `Z_q` 원소이므로 `F` 원소(좌표 0만 채움)로 넣는다. `w_hat`의 자릿수 규약(Greyhound 표의 "table rule": 마지막 자릿수가 더 넓을 수 있음)이 `LevelStatement`의 `G_w`와 다르면 `LevelStatement`를 일반화하되 기존 동작은 바꾸지 마라.
- 제약식을 조밀 벡터로 하나씩 만들지 마라(41개 × 170 MB). 블록 구조를 이용해 집계 벡터 `Phi` 하나만 만든다.
- 공개 행렬 `A, B, D`는 `tail.level_keys`와 같은 방식으로 시드에서 유도한다.
- Fiat–Shamir: `Transcript`에 **공개 파라미터(모양, 랭크, 밑, 버린 비트 수, 모든 레벨의 schedule), 커밋먼트 `u`, `x`, `y`**를 첫 메시지 전에 흡수시키고, first fold의 transcript 상태를 tail이 이어받게 하라. 현재 `tail.prove/verify`는 `bind_statement`로 최상위 statement와 schedule을 흡수한다. PCS에서는 statement가 `(u, x, y, v, y_hat, c)`에서 결정되므로 조밀한 statement를 해시할 필요가 없도록 `prove/verify`에 "이미 묶인 transcript를 받는" 경로를 추가하라.
- 인터페이스: `setup(params)`, `commit(f) -> (u, state)`, `open(state, x) -> (y, proof_bytes)`, `verify(u, x, y, proof_bytes) -> bool`. 증명 바이트에는 `v`, `y_hat`의 3좌표, nonce, 그리고 tail 증명이 들어간다. 바이트 분류(`Writer.sizes`)에 `first_fold` 항목을 추가하라.

파라미터는 모델에서 가져온다.

- 최적화(기대 노름) 설정: `gh.instance(logN, True, 'A')` + `W.search(g['W'], g['beta'], W.WOpt(acc='A'))`, `W.set_bkz(483)`.
- 보장(completeness 증명) 설정: `gh.instance_gc(logN, theta, tmar, lb0, lb)` + `W.WOpt(acc='A', gc=True)`. 이 설정에서는 first fold에서도 `||z||^2 + ||e||^2 <= theta*tau*(B_s^2 + r*n*d*4^(rho_t-1))`가 될 때까지 챌린지를 다시 뽑는다(논문 6.6절 마지막 문단, `g['qmax']`).
- `gh.SHAPE`에는 `N = 2^26, 2^28, 2^30`만 있다. 작은 크기 `2^18, 2^20, 2^22, 2^24`용 모양을 추가하라: `m*r = N/64`, `m/r ≈ 9.5`, `b_0 = 16`(`delta_0 = 8`), 커밋 자릿수 밑 `2^6`. 최적화 설정에서 랭크 고정점이 수렴하지 않으면 원인을 보고하라.
- 모델이 예측하는 first fold 몫(`gh.own_bits`)과 구현 바이트를 비교하라.

## 3. 작업 B — 실험 (`proto/run_e2e.py` 신규, 결과 `proto/e2e.json`)

각 `(N, 설정, 입력 종류)` 조합에 대해 여러 번 실행하고 아래를 기록한다.

| 기록 항목 | 설명 |
|---|---|
| 최종 증명 바이트 | 전체, first fold 몫, 레벨별; 모델 값과의 차이(%) |
| 레벨별 `||w_j|| / beta_j` | first fold 출력(레벨 1 입력) 포함. `tail.STATS` 참고 |
| 레벨별 시도 횟수 | 투영(`jl_tries`), 챌린지(`c_tries`), first fold 챌린지; 평균·최대·히스토그램 |
| 실패 | nonce 소진 등 예외를 잡아 횟수와 레벨을 기록 |
| 시간·메모리 | commit / open / verify 각각의 wall time, 최대 RSS(`resource.getrusage`) |

- 크기: Python으로 `N = 2^18, 2^20, 2^22, 2^24`는 반드시, `2^26`은 가능하면. 각 조합 실행 횟수는 작은 크기 20회, 큰 크기 5회 이상(서로 다른 다항식 시드와 평가점).
- 설정: 증명 기반/483(최적화)과 보장/483 두 가지.
- 입력 종류: (i) `[0, q)` 균등 난수 계수, (ii) 희소(1%만 0이 아님), (iii) 모든 계수가 같은 값, (iv) 짧은 주기로 반복되는 계수, (v) 8비트 작은 계수, (vi) 모든 자릿수가 최댓값인 최악 입력(보장 설정의 한계 `B_s`에 닿는 다항식).
- 음성 테스트: 틀린 `y`, 틀린 `x`, 다른 다항식의 커밋먼트, 증명의 임의 비트 반전 60회, 절단, 덧붙이기 — 모두 거부되어야 한다.
- 기존 tail 실험도 최종 코드로 다시 돌려라: `proto/run_tests.py`, `proto/run_retries.py retries.json 20`, `proto/run_full.py <key> full_X.json`(키는 `README.md` 참고). 가능하면 full-size tail 실행을 최적화 설정에서 20회 반복해 시도 횟수 분포를 구하라.

보고할 때는 "최적화 설정에서 실제 first fold witness의 노름이 모델의 기대 노름 한계 대비 얼마였는가"를 한 문장으로 답할 수 있어야 한다. 한계를 넘는 경우가 있으면 숨기지 말고 빈도와 영향을 적어라.

## 4. 작업 C — 큰 크기를 위한 성능 작업 (선택이지만 권장)

`N = 2^28, 2^30`까지 돌리려면 NumPy schoolbook 곱셈(`ring.py`의 `mul`, `mul_small`, `matvec_small`)으로는 부족하다.

- `matvec_small`(작은 계수 × 균등 행렬)과 `mul`을 C 확장(cffi/ctypes) 또는 numba로 옮겨라. `q = 2^32 - 99`에서는 `X^64 + 1`이 완전 분해되지 않으므로 NTT를 쓰려면 보조 소수들에 대한 CRT가 필요하다. 어떤 방법이든 **기존 `ring.mul`과 비트 단위로 같은 결과**임을 무작위 테스트로 확인하라.
- 메모리: `N = 2^30`이면 다항식만 4 GB다. `s_i`를 int64로 모두 들고 있지 말고 블록 단위로 스트리밍하라(커밋은 블록별, `z`는 `z += c_i * s_i` 누적, 자릿수는 필요할 때 재계산하거나 int8로 저장).
- 실제로 돌린 크기만 보고하라. 검증 시간을 Greyhound·Akita와 비교하는 문장은 같은 기계에서 잰 값이 없으면 쓰지 마라.

## 5. 작업 D — 논문 반영

- `tex/07_impl.tex`: "Scope" 문단을 실제로 한 일에 맞게 고치고, "End-to-end runs" 문단과 표(크기 대 모델, 레벨별 노름·시도 횟수, 입력 종류별 결과, 시간·메모리)를 추가하라. 표 내용은 `build.py`가 `proto/e2e.json`에서 생성한다.
- 초록의 "the first fold is not implemented", 1.2절 Limitations, 8절 "Completeness and real witnesses" 문단을 결과에 맞게 수정하라. 결과가 논문의 주장과 어긋나면(예: 최적화 설정에서 실패가 잦음) 주장을 결과에 맞춰 고쳐라.
- `python3 build.py`로 컴파일해 `Overfull`이 4pt를 넘지 않는지, 참조·인용 경고가 없는지 확인하라.

## 6. 완료 기준

1. `proto/pcs.py`의 commit/open/verify가 위 크기에서 통과하고 음성 테스트가 모두 거부된다.
2. `proto/e2e.json`과 요약 표가 있고, 모델 대비 바이트 차이와 그 원인이 설명되어 있다.
3. 기존 테스트와 재실행 결과가 최종 코드 기준으로 갱신되어 있다.
4. 논문이 새 결과로 다시 빌드되고, 돌리지 못한 항목은 Limitations에 남아 있다.
5. `NOTES.md`에 무엇을 했는지, 논문과 코드 사이에서 발견한 차이, 남은 일이 정리되어 있다.

마지막에 다음을 보고하라: 실제로 돌린 `(N, 설정, 입력)` 목록, 실패 횟수, 모델 대비 최대 편차, 최적화 설정에서 관찰된 `||w||/beta`의 범위와 최대 시도 횟수, 논문에서 고친 문장 목록.

## 7. 알려진 주의점

- `tail.random_instance`는 반올림 때문에 노름이 한계를 1–2% 넘는 witness를 만들 수 있다. 실험에는 `run_full.py`나 `run_retries.py`의 생성기처럼 한계 이하임을 확인하는 생성기를 써라.
- 보장 설정에서 witness 노름은 한계 `B` 이하여야 하고(`tail.STATS`의 `fill <= 1`), 그렇지 않으면 구현이나 모델의 한계식 중 하나가 틀린 것이다. 식은 논문 6.6절의 (13)–(14)(정리 6.3)와 `wmodel._wlevel_gc`, `gh.first_fold_gc`에 있다.
- `model/final.json`의 `trace`가 곧 schedule이다(`run_full.py` 참고). `bjl`은 검증기가 투영에 쓰는 한계, `qmax`는 `(z, e)`의 결합 한계, `joint`는 종단 레벨에서 결합 검사를 쓴다는 표시다.
- 모델 계산(`model/run_final.py`)은 2코어에서 수 시간 걸린다. 새 크기의 schedule은 `W.search`를 직접 불러 구하고 JSON으로 캐시하라.
