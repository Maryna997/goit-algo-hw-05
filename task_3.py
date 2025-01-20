import os
import timeit
from typing import Callable


def read_file(filename):
    project_dir = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(project_dir, filename)
    with open(full_path, "r", encoding="utf-8", errors="ignore") as file:
        return file.read()


def measure_time(algorithm: Callable, text: str, pattern: str):
    setup = f"from __main__ import {algorithm}, read_file; text = read_file('{text}'); pattern = '{pattern}'"
    stmt = f"{algorithm}(text, pattern)"

    time_taken = timeit.timeit(stmt, setup, number=1)
    return time_taken


def boyer_moore(text, pattern):
    m = len(pattern)
    n = len(text)
    if m == 0:
        return 0
    bad_char = {}
    for i in range(m):
        bad_char[pattern[i]] = i
    s = 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
        if j < 0:
            return s
        else:
            s += max(1, j - bad_char.get(text[s + j], -1))
    return -1


def knuth_morris_pratt(text, pattern):
    def compute_lps(pattern):
        lps = [0] * len(pattern)
        length = 0
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    m = len(pattern)
    n = len(text)
    lps = compute_lps(pattern)
    i = 0
    j = 0
    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1
        if j == m:
            return i - j
        elif i < n and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return -1


def rabin_karp(text, pattern, q=101):
    d = 256
    m = len(pattern)
    n = len(text)
    p = 0
    t = 0
    h = 1
    for i in range(m - 1):
        h = (h * d) % q
    for i in range(m):
        p = (d * p + ord(pattern[i])) % q
        t = (d * t + ord(text[i])) % q
    for i in range(n - m + 1):
        if p == t:
            if text[i : i + m] == pattern:
                return i
        if i < n - m:
            t = (d * (t - ord(text[i]) * h) + ord(text[i + m])) % q
            if t < 0:
                t = t + q
    return -1


text_file_1 = "article_1.txt"
text_file_2 = "article_2.txt"
real_pattern_1 = "Це зменшує поле пошуку та робить лінійний пошук життєздатним варіантом"
real_pattern_2 = "Проте для кожного елементу потрібно виділяти новий блок пам’яті."
fake_pattern = "випадковий підрядок"

time_kmp_1 = measure_time("knuth_morris_pratt", text_file_1, real_pattern_1)
time_boyer_moore_1 = measure_time("boyer_moore", text_file_1, real_pattern_1)
time_rabin_karp_1 = measure_time("rabin_karp", text_file_1, real_pattern_1)

dict_real_1 = {
    time_kmp_1: "Knuth-Morris-Pratt",
    time_boyer_moore_1: "Boyer-Moore",
    time_rabin_karp_1: "Rabin-Karp",
}

print()
print(f"Час виконання алгоритмів для існуючого підрядка в {text_file_1}:")
for key, value in dict_real_1.items():
    print(f"{value}: {key}")

time_fake_kmp_1 = measure_time("knuth_morris_pratt", text_file_1, fake_pattern)
time_fake_boyer_moore_1 = measure_time("boyer_moore", text_file_1, fake_pattern)
time_fake_rabin_karp_1 = measure_time("rabin_karp", text_file_1, fake_pattern)

dict_fake_1 = {
    time_fake_kmp_1: "Knuth-Morris-Pratt",
    time_fake_boyer_moore_1: "Boyer-Moore",
    time_fake_rabin_karp_1: "Rabin-Karp",
}

print()
print(f"Час виконання алгоритмів для вигаданого підрядка в {text_file_1}:")
for key, value in dict_fake_1.items():
    print(f"{value}: {key}")

fastest_algorithm_real_1 = min(time_kmp_1, time_boyer_moore_1, time_rabin_karp_1)
fastest_algorithm_fake_1 = min(
    time_fake_kmp_1, time_fake_boyer_moore_1, time_fake_rabin_karp_1
)

dict_best_time = {
    fastest_algorithm_real_1: dict_real_1[fastest_algorithm_real_1]
    + f" (real, {text_file_1})",
    fastest_algorithm_fake_1: dict_fake_1[fastest_algorithm_fake_1]
    + f" (fake, {text_file_1})",
}

print()
print(
    f"Найшвидший алгоритм для існуючого підрядка в {text_file_1}: {dict_real_1[fastest_algorithm_real_1]}: {fastest_algorithm_real_1}"
)
print(
    f"Найшвидший алгоритм для для вигаданого підрядка в {text_file_1}: {dict_fake_1[fastest_algorithm_fake_1]}: {fastest_algorithm_fake_1}"
)

time_kmp_2 = measure_time("knuth_morris_pratt", text_file_2, real_pattern_2)
time_boyer_moore_2 = measure_time("boyer_moore", text_file_2, real_pattern_2)
time_rabin_karp_2 = measure_time("rabin_karp", text_file_2, real_pattern_2)

dict_real_2 = {
    time_kmp_2: "Knuth-Morris-Pratt",
    time_boyer_moore_2: "Boyer-Moore",
    time_rabin_karp_2: "Rabin-Karp",
}

print()
print(f"Час виконання алгоритмів для існуючого підрядка в {text_file_2}:")
for key, value in dict_real_2.items():
    print(f"{value}: {key}")

time_fake_kmp_2 = measure_time("knuth_morris_pratt", text_file_2, fake_pattern)
time_fake_boyer_moore_2 = measure_time("boyer_moore", text_file_2, fake_pattern)
time_fake_rabin_karp_2 = measure_time("rabin_karp", text_file_2, fake_pattern)

dict_fake_2 = {
    time_fake_kmp_2: "Knuth-Morris-Pratt",
    time_fake_boyer_moore_2: "Boyer-Moore",
    time_fake_rabin_karp_2: "Rabin-Karp",
}

print()
print(f"Час виконання алгоритмів для вигаданого підрядка в {text_file_2}:")
for key, value in dict_fake_2.items():
    print(f"{value}: {key}")

fastest_algorithm_real_2 = min(time_kmp_2, time_boyer_moore_2, time_rabin_karp_2)
fastest_algorithm_fake_2 = min(
    time_fake_kmp_2, time_fake_boyer_moore_2, time_fake_rabin_karp_2
)

dict_best_time[fastest_algorithm_real_2] = (
    dict_real_2[fastest_algorithm_real_2] + f" (real, {text_file_2})"
)
dict_best_time[fastest_algorithm_fake_2] = (
    dict_fake_2[fastest_algorithm_fake_2] + f" (fake, {text_file_2})"
)

print()
print(
    f"Найшвидший алгоритм для існуючого підрядка в {text_file_2}: {dict_real_2[fastest_algorithm_real_2]}: {fastest_algorithm_real_2}"
)
print(
    f"Найшвидший алгоритм для для вигаданого підрядка в {text_file_2}: {dict_fake_2[fastest_algorithm_fake_2]}: {fastest_algorithm_fake_2}"
)

min_value_key = min(dict_best_time, key=lambda k: dict_best_time[k])
min_value = dict_best_time[min_value_key]

print()
print(
    f"Найшвидший алгоритм серед усіх замірів виявився алгоритм {min_value}: {min_value_key}"
)