import argparse
import glob
import os
import shutil
from collections import defaultdict
from dataclasses import dataclass
from typing import List
from typing import Dict


def get_file_name_list(all_uniq_dir: str) -> List[str]:
    assert os.path.isdir(all_uniq_dir), f"{all_uniq_dir} is not found!"

    file_name_list: List[str] = [i for i in glob.glob(
        os.path.join(all_uniq_dir, "*.txt"))]
    return file_name_list


def sample_map(map_file: str) -> Dict[str, str]:
    """

    :param map_file:
    :return:  dictionary {R1F1 : MAL-I-1, ... R2F14 : Naive_2}
    """
    assert os.path.isfile(map_file), f"{map_file} is not found!"

    result_dict: Dict[str, str] = dict()
    with open(map_file)as fin:
        fin.readline()
        for line in fin:
            # line looks like this. MAL-I,1,F1,R2
            line: str = line.strip()

            tmp_line_list: List[str] = line.split(",")
            # MAL-I_1
            sample_name: str = f"{tmp_line_list[0]}_{tmp_line_list[1]}"

            # R1F1
            for_rev_name: str = tmp_line_list[3] + tmp_line_list[2]
            # {R1F1 : MAL-I_1}
            result_dict[for_rev_name] = sample_name

    return result_dict


def rename_txt_file(map_dict: Dict[str, str], file_name_list: List[str],
                    rename_dir: str):
    os.makedirs(rename_dir, exist_ok=True)

    for rev_for_name in map_dict:
        for file_name in file_name_list:
            # file_name e.g. 20220329-R1F18-06WIooVL1-2VT17.txt
            file_rf_name: str = os.path.basename(file_name).split("-")[1]
            # get R1F18
            if rev_for_name == file_rf_name:
                print(f"{rev_for_name} matches {file_name}!")
                print(f"rename {file_name} -> {map_dict[rev_for_name]}")
                print()
                shutil.copy(file_name,
                            os.path.join(rename_dir,
                                         map_dict[rev_for_name] + ".txt"))


@dataclass
class SDBInfo:
    ## it looks like this.
    # SDB,Pos,seq,modi_id
    # SDB135,[7:27 52:96],CTTCTGTTTGCTATACCTCTAAGTGTGGAGAAGAATGACCAAAAGACATATCACGCAGGTGGGGGA,1
    sdb_id: str
    pos: str
    seq: str
    modi_id: str
    prefix_pos_1: int
    prefix_pos_2: int
    postfix_pos_1: int
    postfix_pos_2: int


def parse_phage_barcode_map_file(phage_barcode_map_file: str) -> List[SDBInfo]:
    result: List[SDBInfo] = list()
    assert os.path.isfile(
        phage_barcode_map_file), f"{phage_barcode_map_file} is not found!"

    with open(phage_barcode_map_file)as fin:
        fin.readline()
        for line in fin:

            if not line.startswith("SDB"):
                break
            # line looks like this
            # SDB135,[7:27 52:96],CTTCTGTTTGCTATACCTCTAAGTGTGGAGAAGAATGACCAAAAGACATATCACGCAGGTGGGGGA,1
            line: str = line.strip()

            tmp_line_list: List[str] = line.split(",")
            sdb_id: str = tmp_line_list[0]
            pos: str = tmp_line_list[1]
            seq: str = tmp_line_list[2]
            modi_id: str = tmp_line_list[3]

            tmp: str = pos.replace("[", "").replace("]", "")
            prefix_pos: str = tmp.split()[0]

            prefix_pos_1: int = int(prefix_pos.split(":")[0])
            prefix_pos_2: int = int(prefix_pos.split(":")[1])

            postfix_pos: str = tmp.split()[1]
            postfix_pos_1: int = int(postfix_pos.split(":")[0])
            postfix_pos_2: int = int(postfix_pos.split(":")[1])

            tmp_ins: SDBInfo = SDBInfo(sdb_id, pos, seq, modi_id, prefix_pos_1,
                                       prefix_pos_2, postfix_pos_1,
                                       postfix_pos_2)
            result.append(tmp_ins)
    return result


def map_dna_cnt_table_2_sbd_info(phage_sdb_barcode_list: List[SDBInfo],
                                 rename_dir: str, output_dir: str):
    assert os.path.isdir(rename_dir), f"{rename_dir} is not found!"
    log_dir: str = "log_dir"

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    file_name_list: List[str] = [i for i in glob.glob(
        os.path.join(rename_dir, "*.txt"))]

    for cnt_file in file_name_list:
        assert os.path.isfile(cnt_file), f"{cnt_file} is not found!"

        output_samp_file_name: str = os.path.join(output_dir,
                                                  os.path.basename(cnt_file))

        log_samp_file_name: str = os.path.join(log_dir,
                                               os.path.basename(cnt_file))

        with open(cnt_file)as fin, open(output_samp_file_name,
                                        'w')as fout, open(log_samp_file_name,
                                                          'w')as log_out:
            for line in fin:
                line: str = line.strip()
                # line looks like this
                # AAAAAACTACTGTTTGCTATACCGCTGGTGGTACCTTTCTATTCTCACTCTAGTGTGGAGAAGAATGATCAGAAGACTTATCATGCGGGTGGAGGT  KKLLFAIPLVVPFYSHSSVEKNDQKTYHAGGG 5761^M
                tmp_line_list: List[str] = line.split()
                cnt_file_seq: str = tmp_line_list[0]
                cnt_file_cnt: str = tmp_line_list[2]

                for sdb_info in phage_sdb_barcode_list:

                    pre_seq: str = cnt_file_seq[
                                   sdb_info.prefix_pos_1 - 1: sdb_info.prefix_pos_2]
                    post_seq: str = cnt_file_seq[
                                    sdb_info.postfix_pos_1 - 1: sdb_info.postfix_pos_2]
                    cnt_seq: str = pre_seq + post_seq

                    if cnt_seq == sdb_info.seq:
                        out_line: str = f"{sdb_info.sdb_id},{cnt_file_cnt}\n"

                        fout.write(out_line)

                        log_out.write(
                            f"{sdb_info.sdb_id},{pre_seq},{post_seq},{cnt_seq},{cnt_file_cnt}\n")


def aggregation(sdb_dir: str, final_output_dir: str):
    assert os.path.isdir(sdb_dir), f"{sdb_dir} is not found!"

    os.makedirs(final_output_dir, exist_ok=True)

    file_name_list: List[str] = [i for i in glob.glob(
        os.path.join(sdb_dir, "*.txt"))]

    for sdb_file in file_name_list:
        parse_dict: Dict[str, List] = defaultdict(list)

        output_samp_file_name: str = os.path.join(final_output_dir,
                                                  os.path.basename(sdb_file))

        with open(sdb_file)as fin:
            for line in fin:
                line: str = line.strip()

                sdb_id_from_sdb_file: str = line.split(",")[0]
                cnt_from_sdb_file: int = int(line.split(",")[1])

                parse_dict[sdb_id_from_sdb_file].append(cnt_from_sdb_file)

        tmp_dict: Dict[str, int] = {k: sum(v) for (k, v) in parse_dict.items()}
        sorted_parse_dict: List[tuple[str, int]] = sorted(tmp_dict.items(),
                                                          key=lambda kv: kv[1],
                                                          reverse=True)
        with open(output_samp_file_name, 'w')as fout:

            for items in sorted_parse_dict:
                sdb_id_from_dict: str = items[0]
                cnt_from_dict: int = items[1]

                fout.write(f"{sdb_id_from_dict},{cnt_from_dict}\n")


def test(bacode_file: str, all_uniq_dir: str, rename_dir: str,
         phage_barcode_modi_map_file: str, sdb_output_dir: str,
         final_output_dir: str):
    samp_map_dict: Dict[str, str] = sample_map(bacode_file)
    file_name_list: List[str] = get_file_name_list(all_uniq_dir)
    rename_txt_file(samp_map_dict, file_name_list, rename_dir)
    sdb_info_list: List[SDBInfo] = parse_phage_barcode_map_file(
        phage_barcode_modi_map_file)
    map_dna_cnt_table_2_sbd_info(sdb_info_list, rename_dir, sdb_output_dir)
    aggregation(sdb_output_dir, final_output_dir)
    print("done!")


def main(bacode_file: str, all_uniq_dir: str,
         phage_barcode_modi_map_file: str, final_output_dir: str):
    rename_dir: str = 'rename_dir'
    sdb_output_dir: str = "sdb_output_dir"
    samp_map_dict: Dict[str, str] = sample_map(bacode_file)
    file_name_list: List[str] = get_file_name_list(all_uniq_dir)
    rename_txt_file(samp_map_dict, file_name_list, rename_dir)
    sdb_info_list: List[SDBInfo] = parse_phage_barcode_map_file(
        phage_barcode_modi_map_file)
    map_dna_cnt_table_2_sbd_info(sdb_info_list, rename_dir, sdb_output_dir)
    aggregation(sdb_output_dir, final_output_dir)
    print("done!")


if __name__ == '__main__':
    sdb_out_dir: str = "sdb_output_dir"
    final_output_dir: str = "final_output_dir"
    all_uniq_dir: str = "AllUniqueFiles"
    rename_dir: str = "rename_dir"
    barcode_file: str = "AllUniqueFiles/BarcodeSampleMap.csv"
    phage_barcode_modi_map_file: str = "AllUniqueFiles/PhageBarcodeModificationMap.csv"

    ## to avoid argument parser
    # test(barcode_file, all_uniq_dir, rename_dir, phage_barcode_modi_map_file,
    #     sdb_out_dir, final_output_dir)

    parser = argparse.ArgumentParser(
        description="Python script for generating DNA count tables.")

    parser.add_argument("--all_uniq_dir",
                        action='store',
                        type=str,
                        required=True,
                        help="AllUniqueFiles dir"
                        )

    parser.add_argument('--barcode_file',
                        action='store',
                        type=str,
                        required=True,
                        help="BarcodeSampleMap.csv",
                        )

    parser.add_argument('--sdb_map_file',
                        action='store',
                        type=str,
                        required=True,
                        help='PhageBarcodeModificationMap.csv',
                        )
    parser.add_argument('--dna_count_table_output_dir',
                        action='store',
                        type=str,
                        required=True,
                        help='DNA count table output directory',
                        )

    args = parser.parse_args()

    main(args.barcode_file, args.all_uniq_dir, args.sdb_map_file,
         args.final_output_dir)
