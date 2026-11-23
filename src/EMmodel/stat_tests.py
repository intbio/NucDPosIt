from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import norm


class DyadPositionTester:
    """
    Класс для проверки гипотез о положении диады с помощью
    отношения правдоподобия и симуляций Монте-Карло

    H0: позиции неразличимы (можно слить позиции)
    H1: позиции значимо различимы (сливать нельзя)
    """

    def __init__(self, errors, alpha=0.05, nsamples=10000, scan_len=5):
        """
        Параметры:
        - errors: объект ErrorsProbs с распределением ошибок
        - alpha: уровень значимости (по умолчанию 0.05)
        - nsamples: количество симуляций Монте-Карло
        - scan_len: радиус сканирования вокруг dyad_pos
        """
        self.errors = errors
        self.alpha = alpha
        self.nsamples = nsamples
        self.scan_len = scan_len

    def compute_log_likelihood(self, dyad_pos, starts, ends):
        """
        Вычисляет логарифм правдоподобия для заданного положения диады

        Возвращает:
        - log_likelihood: float
        """
        len_left = dyad_pos - starts + 1
        len_right = ends - dyad_pos + 1

        # Защита от выхода за границы
        mask_left = (len_left >= 0) & (len_left < len(self.errors))
        mask_right = (len_right >= 0) & (len_right < len(self.errors))

        logL = 0.0
        if np.any(mask_left):
            logL += np.sum(np.log(self.errors[len_left[mask_left]] + 1e-30))
        if np.any(mask_right):
            logL += np.sum(np.log(self.errors[len_right[mask_right]] + 1e-30))

        return logL

    @lru_cache
    def generate_null_distribution(self, dpos, ntemplates):
        """
        Генерирует распределение статистики критерия при H0
        с помощью симуляций Монте-Карло

        Параметры:
        - dpos: сдвиг между dyad_pos и test_pos
        - ntemplates: количество шаблонов (длина выборки)

        Возвращает:
        - criterion_MK: массив статистик для симуляций
        - H0_MK_sym: массив статистик при H0
        """
        all_templates = ntemplates * self.nsamples

        # Генерация случайных длин
        random_lefts = self.errors.generate_random_len(all_templates).reshape(
            self.nsamples, ntemplates
        )
        random_rights = self.errors.generate_random_len(all_templates).reshape(
            self.nsamples, ntemplates
        )

        # H0: без сдвига
        H0_errors_left = self.errors[random_lefts]
        H0_errors_right = self.errors[random_rights]
        H0_MK_sym = np.sum(np.log(H0_errors_left * H0_errors_right + 1e-30), axis=1)

        # H1: со сдвигом dpos
        H1_random_left = random_lefts - dpos
        H1_left_errors = np.zeros_like(H1_random_left, dtype="float")
        mask = H1_random_left < len(self.errors)
        H1_left_errors[mask] = self.errors[H1_random_left[mask]]

        H1_random_right = random_rights + dpos
        H1_right_errors = np.zeros_like(H1_random_right, dtype="float")
        mask = H1_random_right < len(self.errors)
        H1_right_errors[mask] = self.errors[H1_random_right[mask]]

        H1_MK_sym = np.sum(np.log(H1_right_errors * H1_left_errors + 1e-30), axis=1)

        criterion_MK = H1_MK_sym - H0_MK_sym

        return criterion_MK, H0_MK_sym

    def test_position(self, dyad_pos, test_pos, starts, ends):
        """
        Проверяет гипотезу H0: диада в test_pos против H1: диада в dyad_pos

        Возвращает:
        - results: dict с результатами теста
        """
        dpos = dyad_pos - test_pos

        if dpos == 0:
            return {
                "reject": False,
                "p_value": 1.0,
                "criterion_stat": 0.0,
                "critical_value": 0.0,
                "dpos": 0,
                "message": "Same position",
            }

        # Вычисление наблюдаемой статистики
        h1_logLH = self.compute_log_likelihood(dyad_pos, starts, ends)
        h0_logLH = self.compute_log_likelihood(test_pos, starts, ends)
        criterion_stat = h1_logLH - h0_logLH

        # Генерация нулевого распределения
        ntemplates = len(starts)
        criterion_MK, _ = self.generate_null_distribution(dpos, ntemplates)

        # Критическое значение для правостороннего теста
        right_critical = np.percentile(criterion_MK, 100 * (1 - self.alpha))

        # p-value
        p_value = np.mean(criterion_MK >= criterion_stat)

        # Решение
        reject_h0 = p_value < self.alpha

        return {
            "reject": reject_h0,
            "p_value": p_value,
            "criterion_stat": criterion_stat,
            "critical_value": right_critical,
            "dpos": dpos,
            "ntemplates": ntemplates,
            "distribution": criterion_MK,
            "test_pos": test_pos,
            "dyad_pos": dyad_pos,
        }

    def test_ends_symmerty(self, starts, ends):

        curve1 = np.array(list(occ.values()))  # дискретные значения
        curve2 = np.array(list(template_occ.values()))

        # Находим все возможные уникальные значения
        all_values = np.unique(np.concatenate([curve1, curve2]))

        # Считаем частоты для каждого значения
        freq1 = np.array([np.sum(curve1 == v) for v in all_values])
        freq2 = np.array([np.sum(curve2 == v) for v in all_values])

        # Хи-квадрат тест
        chi2_stat, p_value, dof, expected = stats.chi2_contingency([freq1, freq2])

    def plot_test_result(self, result):
        """
        Визуализирует результат теста с гистограммой и критической областью
        """
        criterion_MK = result["distribution"]
        criterion_stat = result["criterion_stat"]
        right_critical = result["critical_value"]
        p_value = result["p_value"]
        test_pos = result["test_pos"]
        dyad_pos = result["dyad_pos"]
        reject_h0 = result["reject"]

        plt.figure(figsize=(12, 6))

        # Гистограмма распределения статистики при H0
        sns.histplot(
            criterion_MK,
            kde=True,
            bins=50,
            stat="density",
            color="skyblue",
            alpha=0.6,
            label="Distribution of criterion under H0",
        )

        # Получаем текущие пределы по оси Y
        y_min, y_max = plt.gca().get_ylim()

        # Закрашиваем ПРАВУЮ критическую область
        x_vals = np.linspace(right_critical, criterion_MK.max(), 100)
        plt.fill_between(
            x_vals,
            y_min,
            y_max,
            color="red",
            alpha=0.4,
            label=f"Right critical region (α = {self.alpha})\n"
            f"Reject H0 if value ≥ {right_critical:.3f}",
        )

        # Вертикальная линия для наблюдаемого значения
        plt.axvline(
            criterion_stat,
            color="blue",
            linestyle="--",
            linewidth=2.5,
            label=f"Observed criterion = {criterion_stat:.3f}",
        )

        # Критическое значение (порог)
        plt.axvline(
            right_critical,
            color="red",
            linestyle="-",
            linewidth=2,
            label=f"Critical value = {right_critical:.3f}",
        )

        # Оформление графика
        decision = "Reject H0" if reject_h0 else "Do not reject H0"
        color = "green" if not reject_h0 else "red"

        plt.title(
            f"Right-tailed test for dyad position\n"
            f"H0: dyad in test_pos={test_pos} | H1: dyad in dyad_pos={dyad_pos}\n"
            f"p-value = {p_value:.4f} | {decision}",
            fontsize=12,
            fontweight="bold",
            color=color,
        )
        plt.xlabel("Criterion = H1_logLH - H0_logLH", fontsize=11)
        plt.ylabel("Density", fontsize=11)
        plt.legend(loc="upper right")
        plt.grid(alpha=0.3)

        plt.tight_layout()
        plt.show()

    def analyze_group(self, dyad_pos, starts, ends, df=True, plot=False):
        """
        Анализирует всю группу диад: проверяет все test_pos в окрестности

        Параметры:
        - dyad_pos: истинное положение диады
        - group: DataFrame с колонками 'start' и 'stop'

        Возвращает:
        - results: список результатов для каждого test_pos
        """
        test_range = np.arange(dyad_pos - self.scan_len, dyad_pos + self.scan_len + 1)

        all_results = []

        for test_pos in test_range:
            result = self.test_position(dyad_pos, test_pos, starts, ends)
            all_results.append(result)

            # Визуализируем результат (если не одинаковые позиции)
            if result["dpos"] != 0 and plot:
                self.plot_test_result(result)

        if df:
            return pd.DataFrame(all_results)
        return all_results

    def find_confidence_interval(self, dyad_pos, starts, ends):
        """
        Находит доверительный интервал для положения диады
        через инверсию критерия

        Возвращает:
        - ci_lower, ci_upper: границы доверительного интервала
        """

        # Диапазон поиска
        y_min = dyad_pos - self.scan_len
        y_max = dyad_pos + self.scan_len
        y_range = np.arange(y_min, y_max + 1)

        # Для каждого y проверяем, отвергается ли H0
        non_rejected = []

        for test_pos in y_range:
            if test_pos == dyad_pos:
                non_rejected.append(test_pos)
                continue

            result = self.test_position(dyad_pos, test_pos, starts, ends)
            if not result["reject"]:
                non_rejected.append(test_pos)

        # Доверительный интервал
        if len(non_rejected) > 0:
            ci_lower = min(non_rejected)
            ci_upper = max(non_rejected)
        else:
            ci_lower = ci_upper = dyad_pos

        return ci_lower, ci_upper, non_rejected

    def analyze_model(self, model, as_df=True):
        df = model.to_df()
        grouped_df = df.groupby("dyads")
        results = []
        for dyad_pos, group in grouped_df:
            cur_res = self.analyze_group(dyad_pos, group["start"], group["stop"], False)
            results.extend(cur_res)
        if as_df:
            return pd.DataFrame(results)
        return results