
import re

def eq_name(name_1, name_2):
    nl1 = re.split(r'[- ]', name_1.lower())
    nl2 = re.split(r'[- ]', name_2.lower())
    for n1i, n2i in zip(nl1, nl2):
        k = min(len(n1i), len(n2i))
        if n1i[:k] != n2i[:k]:
            return False
    return True 



def normalize_institution(text: str) -> str:
    stop_words = {
        "the", "of", "and", "for", "in", "at", "on", "to", "a", "an"
    }
    if not text:
        return ""
    text = text.lower()
    text = re.sub(r"[^\w\s]", " ", text)
    tokens = [t for t in text.split() if t not in stop_words]
    return " ".join(tokens)

def eq_affiliation(res_inst:str, auth_inst:str) -> bool:
    """This function is used when comparing researcher institution in db and affiliation from outsource"""
    res_insts = re.split(r'[- ,.]', res_inst.lower())
    auth_insts = re.split(r'[- ,.]', auth_inst.lower())
    for inst_part in res_insts.copy():
        if inst_part not in auth_insts:
            return False
    return True


def eq_institution(db_inst:str, new_inst:str) -> bool:
    """This function is used when comparing institution in db and new institution"""
    norm_db_inst = normalize_institution(db_inst)
    norm_new_inst = normalize_institution(new_inst)
    return norm_db_inst == norm_new_inst


