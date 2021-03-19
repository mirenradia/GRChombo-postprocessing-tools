#!/usr/bin/env bash

function print_help {
  echo "Usage: $0 [OPTION]... [INPUT FILES]..."
  echo "Options:"
  echo "-d NUM_COLUMNS    drop the first NUM_COLUMNS columns from each file"
  echo "                  after the first"
  echo "                  default: 1"
  echo "-s START_WIDTH    width of the first NUM_COLUMNS column in the file"
  echo "                  default: 13"
  echo "-w WIDTH          widths of the other column in the output"
  echo "                  default: 15"
  echo "-c                include comment lines"
  echo "                  default: false"
  echo "-h                print help"
}

OPTSTRING="d:s:w:ch"
# Default number of columns to drop
DROP_COLUMNS=1
# Default to dropping comment lines
INCLUDE_COMMENT_LINES=false
# Default to a column width of 13 for the first DROP_COLUMNS columns
START_COLUMN_WIDTH="13"
# Default to a column width of 15 for the other columns
COLUMN_WIDTH="15"

while getopts ${OPTSTRING} ARG;
do
  case ${ARG} in
    d)
      if ! [[ ${OPTARG} =~ ^[0-9]+$ ]]; then
        echo "argument to -d must be a non-negative integer"
        print_help
        exit 1
      fi
      DROP_COLUMNS=${OPTARG}
      ;;
    s)
      if ! [[ ${OPTARG} =~ ^[0-9]+$ ]]; then
        echo "argument to -w must be a non-negative integer"
        print_help
        exit 1
      fi
      START_COLUMN_WIDTH=${OPTARG}
      ;;
    w)
      if ! [[ ${OPTARG} =~ ^[0-9]+$ ]]; then
        echo "argument to -w must be a non-negative integer"
        print_help
        exit 1
      fi
      COLUMN_WIDTH=${OPTARG}
      ;;
    c)
      INCLUDE_COMMENT_LINES=true
      ;;
    \?)
      print_help
      exit 1
      ;;
    h)
      print_help
      exit 0
      ;;
  esac
done

# Get all the files to join
ALLARGS=("$@")
FILES=("${ALLARGS[@]:($OPTIND-1)}")
echo -n "# Joined from: "
for FILE in "${FILES[@]}";
do
  echo -n "$FILE "
done
# Print new line
echo ""

# Check at least 1 file has been provided
if [[ ${#FILES[@]} -lt 1 ]]; then
  echo "Please provide at least 1 file";
  exit 2;
fi

# Check files exist
for FILE in "${FILES[@]}"
do
  if ! [[ -f "$FILE" ]]; then
    echo "$FILE does not exist."
    exit 3
  fi
done

# Check all files have the same length
if $INCLUDE_COMMENT_LINES; then
  NUM_LINES=$(wc -l ${FILES[0]})
else
  NUM_LINES=$(sed '/^\s*#.*$/ d' "${FILES[0]}" | wc -l)
fi
# Make some temporary files without the comment lines
NOCOMMENT_FILES=()
for FILE in "${FILES[@]}"
do
  if ! $INCLUDE_COMMENT_LINES; then
    NOCOMMENT_FILE="/tmp/nocomment-$(basename ${FILE})"
    NOCOMMENT_FILES+=("$NOCOMMENT_FILE")
    sed '/^\s*#.*$/ d' $FILE > "$NOCOMMENT_FILE"
    FILE_TO_USE=$NOCOMMENT_FILE
  else
    FILE_TO_USE=$FILE
  fi
  NUM_LINES_CHECK=$(wc -l "$FILE_TO_USE" | awk '{print $1}')
  if [[ "$NUM_LINES_CHECK" != "$NUM_LINES" ]]; then
    echo "$FILE has a different number of lines to ${FILES[0]}"
    exit 4
  fi
done

if $INCLUDE_COMMENT_LINES; then
  FILES_TO_USE=("${FILES[@]}")
else
  FILES_TO_USE=("${NOCOMMENT_FILES[@]}")
fi

# Remove DROP_COLUMNS columns from all but the first file and write these
# to temporary files
DECOLUMNED_FILES_TO_USE=("${FILES_TO_USE[0]}")
for FILE in "${FILES_TO_USE[@]:1}"
do
  DECOLUMNED_FILE="/tmp/decolumned-$(basename ${FILE})"
  awk -v drop_columns=${DROP_COLUMNS} \
  '{
    for (i = 1; i <= drop_columns; ++i) {
      $i = ""
   }; print}' $FILE > "$DECOLUMNED_FILE"
   DECOLUMNED_FILES_TO_USE+=("${DECOLUMNED_FILE}")
done

# Now paste the columns together and print nicely with awk
paste "${DECOLUMNED_FILES_TO_USE[@]}" | awk \
  -v start_column_width=${START_COLUMN_WIDTH} -v column_width=${COLUMN_WIDTH} \
  -v drop_columns=${DROP_COLUMNS} \
    '{
      for (i = 1; i <= drop_columns; ++i){
        printf "%*s", start_column_width, $i}
      for (j = drop_columns + 1; j <= NF; ++j){
      printf "%*s", column_width, $j}
      printf "\n"
    }'

# Now clean up temporary files
rm -f "${NOCOMMENT_FILE[@]}"
rm -f "${DECOLUMNED_FILES_TO_USE[@]}"
