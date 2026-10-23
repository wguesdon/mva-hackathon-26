-- A two-column layout rejects the longtable that Quarto writes for a markdown table. This filter
-- writes each table as a booktabs tabular instead, and it wraps the text of a wide column so that
-- the table fits the measure.
--
-- Three cases follow from the width of the widest cell of each column, measured in characters.
--
--   1. A table that fits the column is a plain tabular.
--   2. A table that is too wide for the column becomes a tabularx of the column width, where every
--      column wider than WRAP_AT characters becomes a wrapping L column and the rest keep their
--      natural width.
--   3. A table that is too wide for the page keeps the same form. The float that carries it is
--      widened to both columns by the report source, which writes the table inside a table*
--      environment.
--
-- The capacities below are characters of the body font at 9 point, measured from the rendered
-- output of this layout.
local COLUMN_CAPACITY = 46
local WRAP_AT = 18

local function cell_latex(cell)
  local doc = pandoc.Pandoc(cell.contents)
  local text = pandoc.write(doc, "latex")
  return (text:gsub("\n+$", ""):gsub("\n", " "))
end

local function all_rows(tbl)
  local rows = {}
  for _, row in ipairs(tbl.head.rows) do
    table.insert(rows, row)
  end
  for _, body in ipairs(tbl.bodies) do
    for _, row in ipairs(body.body) do
      table.insert(rows, row)
    end
  end
  return rows
end

-- The longest rendered cell of each column, measured in characters.
local function column_widths(tbl)
  local widths = {}
  for index = 1, #tbl.colspecs do
    widths[index] = 0
  end
  for _, row in ipairs(all_rows(tbl)) do
    for index, cell in ipairs(row.cells) do
      local length = #cell_latex(cell)
      if widths[index] and length > widths[index] then
        widths[index] = length
      end
    end
  end
  return widths
end

-- Three characters of padding sit beside every column.
local function natural_width(widths)
  local total = 0
  for _, width in ipairs(widths) do
    total = total + width + 3
  end
  return total
end

local function align_letter(colspec)
  local map = { AlignLeft = "l", AlignRight = "r", AlignCenter = "c", AlignDefault = "l" }
  return map[colspec[1]] or "l"
end

-- The column specification of the table, and whether the table needs tabularx.
local function build_spec(tbl, widths, wide)
  local spec = {}
  local wrapping = false
  for index, colspec in ipairs(tbl.colspecs) do
    if wide and widths[index] and widths[index] > WRAP_AT then
      table.insert(spec, "L")
      wrapping = true
    else
      table.insert(spec, align_letter(colspec))
    end
  end
  return table.concat(spec), wrapping
end

-- Permit a line break at existing punctuation without changing identifier characters.
-- The preamble already permits a break after an underscore.
local function breakable_code(group)
  return "\\texttt" .. group:gsub("([/%.%-])", "%1\\allowbreak{}")
end

local function breakable_latex(text)
  return (text:gsub("\\texttt(%b{})", breakable_code))
end

local function row_latex(row)
  local cells = {}
  for _, cell in ipairs(row.cells) do
    table.insert(cells, breakable_latex(cell_latex(cell)))
  end
  return table.concat(cells, " & ") .. " \\\\"
end

local function header_latex(row)
  local cells = {}
  for _, cell in ipairs(row.cells) do
    table.insert(cells, "\\textsf{\\bfseries " .. cell_latex(cell) .. "}")
  end
  return table.concat(cells, " & ") .. " \\\\"
end

function Table(tbl)
  if not FORMAT:match("latex") then
    return nil
  end
  local widths = column_widths(tbl)
  local wide = natural_width(widths) > COLUMN_CAPACITY
  local spec, wrapping = build_spec(tbl, widths, wide)
  local lines = {}
  table.insert(lines, "\\begingroup\\footnotesize")
  if wrapping then
    table.insert(lines, "\\begin{tabularx}{\\linewidth}{" .. spec .. "}")
  else
    table.insert(lines, "\\begin{tabular}{" .. spec .. "}")
  end
  table.insert(lines, "\\toprule")
  for _, row in ipairs(tbl.head.rows) do
    table.insert(lines, header_latex(row))
  end
  table.insert(lines, "\\midrule")
  for _, body in ipairs(tbl.bodies) do
    for _, row in ipairs(body.body) do
      table.insert(lines, row_latex(row))
    end
  end
  table.insert(lines, "\\bottomrule")
  if wrapping then
    table.insert(lines, "\\end{tabularx}")
  else
    table.insert(lines, "\\end{tabular}")
  end
  table.insert(lines, "\\endgroup")
  return pandoc.RawBlock("latex", table.concat(lines, "\n"))
end
