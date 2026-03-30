import argparse
import pandas as pd
import os

# Definição das chaves compostas que identificam uma simulação única
CHAVES_UNICAS = ['fibrosis_type', 'fiber_angle_deg', 'density', 'seed', 'atpi']

def sort_dataframe(df):
    """
    Ordena o DataFrame:
    1. Ordem alfabética de fibrosis_type
    2. Ordem crescente de fiber_angle_deg
    3. Ordem crescente de density
    """
    # Garantir os tipos de dados para ordenação correta
    df['fibrosis_type'] = df['fibrosis_type'].astype(str)
    df['fiber_angle_deg'] = pd.to_numeric(df['fiber_angle_deg'])
    df['density'] = pd.to_numeric(df['density'])

    df_sorted = df.sort_values(
        by=['fibrosis_type', 'fiber_angle_deg', 'density'],
        ascending=[True, True, True]
    )
    return df_sorted.reset_index(drop=True)

def rename_stochastic_to_uniform(df):
    """Substitui ocorrências de 'stochastic' (ou 'Stochastic') por 'uniform/Uniform'."""
    # Trata tanto minúsculas quanto a primeira maiúscula
    df['fibrosis_type'] = df['fibrosis_type'].replace('stochastic', 'uniform')
    df['fibrosis_type'] = df['fibrosis_type'].replace('Stochastic', 'Uniform')
    print("[Ação] Todas as ocorrências de 'stochastic' foram alteradas para 'uniform'.")
    return df

def fill_missing(df_1, df_2):
    """
    Acrescenta no Arquivo 1 apenas as chaves que existem no Arquivo 2, mas não no Arquivo 1.
    As linhas originais do Arquivo 1 permanecem inalteradas.
    """
    df1_idx = df_1.set_index(CHAVES_UNICAS)
    df2_idx = df_2.set_index(CHAVES_UNICAS)

    # Identifica índices que estão no 2 mas não estão no 1
    missing_indices = df2_idx.index.difference(df1_idx.index)
    df_missing = df2_idx.loc[missing_indices]

    # Concatena as linhas novas ao df1
    df_result = pd.concat([df1_idx, df_missing]).reset_index()
    print(f"[Ação] Preenchimento concluído. {len(df_missing)} novas simulações foram adicionadas do Arquivo 2 para o Arquivo 1.")
    return df_result

def interactive_merge(df_1, df_2):
    """
    Adiciona linhas do Arquivo 2 ao Arquivo 1.
    Se a chave existir em ambos com 'final_time_ms' diferente, pergunta ao usuário.
    """
    df1_idx = df_1.set_index(CHAVES_UNICAS)
    df2_idx = df_2.set_index(CHAVES_UNICAS)

    # 1. Identificar chaves exclusivas do arquivo 2 (adicioná-las diretamente)
    missing_indices = df2_idx.index.difference(df1_idx.index)

    # 2. Identificar chaves comuns para verificar conflitos
    common_keys = df1_idx.index.intersection(df2_idx.index)

    conflitos_resolvidos = 0

    for key in common_keys:
        val_1 = df1_idx.loc[key, 'final_time_ms']
        val_2 = df2_idx.loc[key, 'final_time_ms']

        # Lidar com duplicatas no mesmo arquivo (segurança)
        if isinstance(val_1, pd.Series): val_1 = val_1.iloc[0]
        if isinstance(val_2, pd.Series): val_2 = val_2.iloc[0]

        if val_1 != val_2:
            print(f"\n[!] CONFLITO ENCONTRADO na chave: {key}")
            print(f"  [1] Valor no Arquivo 1: {val_1}")
            print(f"  [2] Valor no Arquivo 2: {val_2}")

            while True:
                escolha = input("-> Qual valor você deseja manter na versão final? (Digite 1 ou 2): ").strip()
                if escolha == '1':
                    break # Mantém o que já está no df1_idx
                elif escolha == '2':
                    df1_idx.loc[key, 'final_time_ms'] = val_2 # Atualiza com o valor do arquivo 2
                    conflitos_resolvidos += 1
                    break
                else:
                    print("Opção inválida. Digite 1 ou 2.")

    # Adicionar as exclusivas do arquivo 2
    df_result = pd.concat([df1_idx, df2_idx.loc[missing_indices]]).reset_index()

    print(f"\n[Ação] Merge interativo concluído. {len(missing_indices)} linhas novas adicionadas.")
    print(f"[Ação] {conflitos_resolvidos} conflitos foram resolvidos substituindo pelos valores do Arquivo 2.")
    return df_result

def main():
    parser = argparse.ArgumentParser(description="Ferramenta modular de manipulação de CSVs de Simulação.")
    parser.add_argument("--file1", required=True, help="Caminho para o CSV principal (Arquivo 1).")
    parser.add_argument("--file2", required=False, help="Caminho para o CSV secundário (Arquivo 2), necessário para --fill e --merge.")
    parser.add_argument("--out", required=False, help="Caminho para o arquivo de saída. Se não fornecido, sobrescreve o file1.")

    # Ações modulares
    parser.add_argument("--rename", action="store_true", help="Substitui 'stochastic' por 'uniform' no Arquivo 1.")
    parser.add_argument("--fill", action="store_true", help="Completa o Arquivo 1 com linhas faltantes presentes no Arquivo 2.")
    parser.add_argument("--merge", action="store_true", help="Funde o Arquivo 2 no 1, perguntando ao usuário em caso de conflitos.")

    args = parser.parse_args()

    # Validações iniciais
    if not os.path.exists(args.file1):
        print(f"Erro: O arquivo {args.file1} não foi encontrado.")
        return

    if (args.fill or args.merge) and not args.file2:
        print("Erro: Para usar --fill ou --merge, você deve fornecer o Arquivo 2 através de --file2.")
        return

    if args.fill and args.merge:
        print("Aviso: Você ativou --fill e --merge. O --merge engloba o --fill (com interatividade). Executando apenas o --merge.")
        args.fill = False

    # Carrega o arquivo principal
    df_main = pd.read_csv(args.file1)

    # Execução sequencial das funções modulares
    if args.fill:
        df_sec = pd.read_csv(args.file2)
        df_main = fill_missing(df_main, df_sec)

    if args.merge:
        df_sec = pd.read_csv(args.file2)
        df_main = interactive_merge(df_main, df_sec)

    if args.rename:
        df_main = rename_stochastic_to_uniform(df_main)

    # Sempre ordenar no final, independentemente das ações
    df_main = sort_dataframe(df_main)

    # Salvar o resultado
    arquivo_saida = args.out if args.out else args.file1
    df_main.to_csv(arquivo_saida, index=False)
    print(f"\nProcessamento finalizado! Arquivo salvo em: {arquivo_saida}")

if __name__ == "__main__":
    main()
