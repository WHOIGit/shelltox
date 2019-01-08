
const blob = new Blob([file_content], { type: 'text/csv;charset=utf-8;' })

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, file_name)
} else {
    const link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = file_name
    link.target = '_blank'
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}
