from mypy_boto3_s3.service_resource import ObjectSummary
import boto3
import json

import logging

logging.getLogger('boto3').setLevel(logging.WARNING)
logging.getLogger('botocore').setLevel(logging.WARNING)

def get_basic_object_metadata(obj: ObjectSummary) -> dict:
    """
    This function retrieves basic metadata of an AWS S3 object.

    Parameters:
        obj (ObjectSummary): The AWS S3 object.

    Returns:
        dict: A dictionary with the size, ETag, last modified date, storage platform, region, and storage tier of the object.
    """
    _ = obj.load()
    return {
        "file:size": obj.content_length,
        "e_tag": obj.e_tag.strip('"'),
        "last_modified": obj.last_modified.isoformat(),
        "storage:platform": "AWS",
        "storage:region": obj.meta.client.meta.region_name,
        "storage:tier": obj.storage_class,
    }


def copy_item_to_s3(item, s3_key):
    """
    This function copies an item to an AWS S3 bucket.

    Parameters:
        item: The item to copy. It must have a `to_dict` method that returns a dictionary representation of it.
        s3_key (str): The file path in the S3 bucket to copy the item to.

    The function performs the following steps:
        1. Initializes a boto3 S3 client and splits the s3_key into the bucket name and the key.
        2. Converts the item to a dictionary, serializes it to a JSON string, and encodes it to bytes.
        3. Puts the encoded JSON string to the specified file path in the S3 bucket.
    """
    s3 = boto3.client("s3")
    bucket, key = split_s3_key(s3_key)
    s3.put_object(Body=json.dumps(item.to_dict()).encode(), Bucket=bucket, Key=key)


def split_s3_key(s3_key: str) -> tuple[str, str]:
    """
    This function splits an S3 key into the bucket name and the key.

    Parameters:
        s3_key (str): The S3 key to split. It should be in the format 's3://bucket/key'.

    Returns:
        tuple: A tuple containing the bucket name and the key. If the S3 key does not contain a key, the second element
          of the tuple will be None.

    The function performs the following steps:
        1. Removes the 's3://' prefix from the S3 key.
        2. Splits the remaining string on the first '/' character.
        3. Returns the first part as the bucket name and the second part as the key. If there is no '/', the key will
          be None.
    """
    parts = s3_key.replace("s3://", "").split("/", 1)
    bucket = parts[0]
    key = parts[1] if len(parts) > 1 else None
    return bucket, key


def s3_key_public_url_converter(url: str) -> str:
    """
    This function converts an S3 URL to an HTTPS URL and vice versa.

    Parameters:
        url (str): The URL to convert. It should be in the format 's3://bucket/' or 'https://bucket.s3.amazonaws.com/'.

    Returns:
        str: The converted URL. If the input URL is an S3 URL, the function returns an HTTPS URL. If the input URL is
        an HTTPS URL, the function returns an S3 URL.

    The function performs the following steps:
        1. Checks if the input URL is an S3 URL or an HTTPS URL.
        2. If the input URL is an S3 URL, it converts it to an HTTPS URL.
        3. If the input URL is an HTTPS URL, it converts it to an S3 URL.
    """
    if url.startswith("s3://"):
        bucket = url.replace("s3://", "").split("/")[0]
        key = url.replace(f"s3://{bucket}", "")[1:]
        return f"https://{bucket}.s3.amazonaws.com/{key}"
    elif url.startswith("https://") and ".s3.amazonaws.com/" in url:
        bucket = url.replace("https://", "").split(".s3.amazonaws.com")[0]
        key = url.replace(f"https://{bucket}.s3.amazonaws.com/","")
        return f"s3://{bucket}/{key}"
    else:
        raise ValueError("Invalid URL format")
    

def verify_safe_prefix(s3_key: str):
    """
    TODO: discuss this with the team. Would like some safety mechanism to ensure that the S3 key is limited to  
    certain prefixes. Should there be some restriction where these files can be written?
    """
    parts = s3_key.split("/")
    if parts[3] != "stac":
        raise ValueError("prefix must begin with stac/")
