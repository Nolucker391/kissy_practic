import math

from sympy import Matrix


def write_string():
    """
    Записывает разделительную строку в файл.
    :return:
    """
    file.write('#' + '-' * 80 + '#\n')


def function_cramer(coefficients, constants):
    """
    Функция крамера для вычисления системы уравнений.

    :coefficients = list(Matrix())
    :constants = list(Matrix())
    """
    coeff_det = coefficients.det()
    solution_list = []

    for index in range(coefficients.shape[0]):
        Matrix_coefficients_index = coefficients.copy()
        Matrix_coefficients_index[:, index] = constants

        solution_list.append(Matrix_coefficients_index.det() / coeff_det)

    return solution_list


with open('result.txt', 'w') as file:
    o = 0
    P = 250
    T = 2700
    R0 = 1.987
    k = 0.8
    b = [2, k]
    p_ = [P for i in range(6)]
    psum = sum(p_)
    Mt = P


    write_string()

    file.write(
        'Исходные данные по варианту №21\n'
        'студента группы: 23СЖ02 Губайдуллина Арслана\n'
    )

    write_string()

    file.write(
        f'Температура[T] = {T}\n'
        f'Давление[P] = {P} кг/см2\n'
        f'Коэффициент избытка окислителя[alpha] = {k}\n'
    )

    I_values = [
        18579.6,
        20605.7,
        -31302.2,
        28228.5,
        64033.6,
        71581
    ]
    S_values = [
        47.5422,
        66.9743,
        67.1675,
        60.4515,
        38.3377,
        49.5693
    ]
    Cp_values = [
        8.71414,
        9.37875,
        13.3858,
        8.74102,
        4.96821,
        4.98706
    ]

    ko = ['H2', 'O2', 'H2O', 'OH', 'H', 'O']
    for i in range(6):
        file.write(f'I_{ko[i]} = {I_values[i]} ккал/кмоль\n')
    for i in range(6):
        file.write(f'S_{ko[i]} = {S_values[i]} ккал/(кмоль*K)\n')
    for i in range(6):
        file.write(f'Cp_{ko[i]} = {Cp_values[i]} ккал/(кмоль*K)\n')

    write_string()

    print('Результаты расчета:', file=file)

    write_string()

    indexs_ = [[2, 0, 2, 1, 1, 0],
               [0, 2, 1, 1, 0, 1]]

    k_k = [0 for i in range(4)]
    for j in range(4):
        sum1 = 0
        sum2 = 0
        for i in range(2):
            sum1 += indexs_[i][j] * S_values[i + 4]
            sum2 += indexs_[i][j] * I_values[i + 4]
            k_k[j] = ((sum1 - S_values[j]) / R0 - (sum2 - I_values[j]) / R0 / T) / 2.31

    k_kk = [0 for i in range(4)]
    for j in range(4):
        sum3 = 0
        for i in range(2):
            sum3 += indexs_[i][j] * I_values[i + 4]
            k_kk[j] = ((sum3 - I_values[j]) / R0 / T)

    while True:
        p_sum = 0
        for i in range(6):
            p_sum += p_[i]

        B_ = [0 for i in range(2)]
        for i in range(2):
            for j in range(6):
                B_[i] += indexs_[i][j] * p_[j]

        s_k = [0 for i in range(6)]
        for j in range(4):
            sum4 = 0
            for i in range(2):
                sum4 += indexs_[i][j] * math.log10(p_[i + 4])
                s_k[j] = math.log10(p_[j]) - sum4 + k_k[j]
                s_k[i + 4] = math.log10(B_[i]) - math.log10(Mt) - math.log10(b[i])
        s_kp = math.log10(p_sum) - math.log10(P)
        coefficients = [[1, 0, 0, 0, -2, 0, 0],
                        [0, 1, 0, 0, 0, -2, 0],
                        [0, 0, 1, 0, -2, -1, 0],
                        [0, 0, 0, 1, -1, -1, 0],
                        [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
                        [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
                        [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]]

        constants = Matrix([-s_k[0], -s_k[1], -s_k[2], -s_k[3], -s_k[4] * (B_[0]), -s_k[5] * (B_[1]), -s_kp * p_sum])

        solution = function_cramer(Matrix(coefficients), constants)

        for i in range(6):
            p_[i] = 10 ** (math.log10(p_[i]) + solution[i])
        Mt = 10 ** (math.log10(Mt) + solution[6])

        svch = sum(constants) - (-s_kp * p_sum)
        o += 1
        print('приближение', o, file=file)

        for i in range(len(coefficients)):
            for j in range(len(coefficients[i])):
                print("{:22}".format(coefficients[i][j]), end=" ", file=file)
            if i < len(solution):
                print("{:.10f}".format(constants[i]), file=file)
        if abs(svch) < 0.000001:
            break


    coefficients_second = [[1, 0, 0, 0, -2, 0, 0],
                           [0, 1, 0, 0, 0, -2, 0],
                           [0, 0, 1, 0, -2, -1, 0],
                           [0, 0, 0, 1, -1, -1, 0],
                           [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
                           [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
                           [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]]
    constants_second = Matrix([-k_kk[0], -k_kk[1], -k_kk[2], -k_kk[3], 0, 0, 0])
    solution_second = function_cramer(Matrix(coefficients_second), constants_second)

    coefficients_third = [[1, 0, 0, 0, -2, 0, 0],
                          [0, 1, 0, 0, 0, -2, 0],
                          [0, 0, 1, 0, -2, -1, 0],
                          [0, 0, 0, 1, -1, -1, 0],
                          [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
                          [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
                          [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]]
    constants_third = Matrix([0, 0, 0, 0, 0, 0, p_sum])
    solution_third = function_cramer(Matrix(coefficients_third), constants_third)

    mut = 1.008 * 2 + k * 16
    sum5 = 0

    for j in range(6):
        sum5 += I_values[j] * p_[j]
    I = 4184 * sum5 / Mt / mut
    sum6 = 0

    for j in range(6):
        sum6 += S_values[j] * p_[j]
    S = 4184 * sum6 / Mt / mut
    mu = Mt * mut / P
    r = 8314 / mu
    alpha = (1 / T) * (1 - solution_second[6])
    beta = (1 / P) * (solution_third[6])
    sum7 = 0

    for j in range(6):
        sum7 += Cp_values[j] * p_[j]
    cpf = 4184 * sum7 / Mt / mut
    kf = 1 / (1 - (r / cpf))
    af = (T * r * kf) ** 0.5
    sum8 = 0

    for j in range(6):
        sum8 += I_values[j] * p_[j] * solution_second[j]
    cpe = cpf + (1 / T) * (4184 * sum8 / Mt / mut - I * solution_second[6])
    kr = cpe / (cpe - ((r * (1 - solution_second[6]) ** 2) / solution_third[6]))
    a = (T * r * kr / solution_third[6]) ** 0.5

    for i in range(6):
        file.write(f'p_{ko[i]} = {p_[i]} кг/см2\n')
    file.write(f'\t Число итераций - {o}\n')
    file.write(
        f'\t "I" = {I} Дж/кг/К\n'
        f'\t "S" = {S} Дж/кг/К\n'
        f'\t "Mu" = {mu} кг/кмоль\n'
        f'\t "R" = {r} Дж/кг/К\n'
        f'\t "alpha_p" = {alpha} 1/К\n'
        f'\t beta_t = {beta} м2/H\n'
        f'\t cp_f = {cpf} Дж/кг/К\n'
        f'\t k_f = {kf}\n'
        f'\t a_f = {af} м/с\n'
        f'\t c_p_e = {cpe} Дж/кг/К\n'
        f'\t k_R = {kr}\n'
        f'\t a = {a} м/с\n'
    )
    write_string()
